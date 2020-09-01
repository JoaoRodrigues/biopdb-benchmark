#!/usr/bin/env python

"""
Script to process PDBML files.
"""

import argparse
import gzip
import os
import pathlib
import re
import sys
import time
import traceback
from xml.dom import minidom  # for pretty-printing
from xml.etree import ElementTree


def main():

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        'folder',
        type=pathlib.Path,
        help='Top-level folder with PDBML files'
    )
    args = ap.parse_args()

    flist = sorted(args.folder.rglob('*.xml.gz'))
    print(f'Found {len(flist)} files')

    for idx, fpath in enumerate(flist, start=1):
        print(f'Parsing file: {fpath.parent.name}/{fpath.name}')

        # Parse
        with gzip.open(fpath, mode='rt') as f:
            tree = ElementTree.parse(f)
            root = tree.getroot()
            try:
                prefix = re.match('(\{.*\})', root.tag)
                prefix = prefix.group(1)
            except Exception as err:
                raise Exception('Could not parse XML namespace') from err

            models = set()
            chains = set()
            resids = set()
            aatoms = set()
            hetatoms = set()
            for atom in tree.findall(f'.//{prefix}atom_site'):

                try:
                    model = atom.find(f'{prefix}pdbx_PDB_model_num').text
                except AttributeError:
                    model = None

                try:
                    chain = atom.find(f'{prefix}auth_asym_id').text
                except AttributeError:
                    chain = atom.find(f'{prefix}label_asym_id').text

                try:
                    resid = atom.find(f'{prefix}auth_seq_id').text
                except AttributeError:
                    resid = atom.find(f'{prefix}label_seq_id').text

                try:
                    resn = atom.find(f'{prefix}auth_comp_id').text
                except AttributeError:
                    resn = atom.find(f'{prefix}label_comp_id').text

                try:
                    icode = atom.find(f'{prefix}pdbx_PDB_ins_code').text
                except AttributeError:
                    icode = ' '

                try:
                    hetatm = atom.find(f'{prefix}group_PDB').text == 'HETATM'
                except AttributeError:
                    hetatm = False  # just in case

                try:
                    altloc = atom.find(f'{prefix}label_alt_id').text
                except AttributeError:
                    altloc = ' '

                try:
                    atname = atom.find(f'{prefix}auth_atom_id').text
                except AttributeError:
                    atname = atom.find(f'{prefix}label_atom_id').text

                models.add(model)
                chains.add(chain)
                resids.add((chain, resid, resn, icode))
                aatoms.add((chain, resid, resn, icode, atname, altloc, hetatm))

            n_models = len(models)
            n_chains = len(chains)
            n_resids = len(resids)
            n_aatoms = len(aatoms)

            n_icodes = 0
            n_ptmuts = 0
            ptmuts = set()
            for r in resids:
                n_icodes += int(r[3] != ' ')
                rid = (r[0], r[1], r[3])
                if rid not in ptmuts:
                    ptmuts.add(rid)
                else:
                    n_ptmuts += 1

            n_hetatm = 0
            n_altloc = 0
            for a in aatoms:
                n_altloc += int(a[-2] != ' ')
                n_hetatm += int(a[-1])

            data = {
                'models': n_models,
                'chains': n_chains,
                'residues': n_resids,
                'res-icode': n_icodes,
                'res-ptmut': n_ptmuts,
                'all-atoms': n_aatoms,
                'het-atoms': n_hetatm,
                'altlocs': n_altloc
            }

            # Write XML file with numbers
            root = ElementTree.Element('structure')
            root.set('path', fpath.name)

            for key, value in data.items():
                child = ElementTree.SubElement(root, key)
                child.text = str(value)

            # Reparse for pretty print
            xml = minidom.parseString(
                ElementTree.tostring(root, 'utf-8')
            )

            # Write to file
            with fpath.with_suffix('.parsed.xml').open('w') as f:
                f.write(xml.toprettyxml(indent='  '))


            root.clear()
            del tree, root, xml


if __name__ == '__main__':
    main()
