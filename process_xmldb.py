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
            n_atoms = 0  # easier
            for atom in tree.findall(f'.//{prefix}atom_site'):
                n_atoms += 1

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
                    icode = atom.find(f'{prefix}pdbx_PDB_ins_code').text
                except AttributeError:
                    icode = None

                models.add(model)
                chains.add(chain)
                resids.add((chain, resid, icode))

            data = {
                'models': len(models),
                'chains': len(chains),
                'resids': len(resids),
                'atoms': n_atoms
            }

            # Write XML file with numbers
            root = ElementTree.Element('structure')
            root.set('path', fpath.name)

            for level in ('models', 'chains', 'resids', 'atoms'):
                child = ElementTree.SubElement(root, level)
                child.text = str(data.get(level))

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
