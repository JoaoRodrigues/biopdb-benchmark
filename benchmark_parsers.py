#!/usr/bin/env python

"""
Script to benchmark Bio.PDB parsers on a collection of files.

Parses all files in a folder (expected gziped by the way),
writes them out, and parses them again.

Outputs XML files with some numbers on the parsed Structure
objects for later comparison with other datasets/references.
"""

import argparse
import gzip
import logging
import os
import pathlib
import sys
import time
import traceback
from xml.dom import minidom  # for pretty-printing
from xml.etree.ElementTree import (
    Element, SubElement, tostring,
)

import psutil

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBIO, MMCIFIO


def setup_logging():

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        fmt='[%(asctime)s] %(message)s',
        datefmt='%H:%M:%S'
    )
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.handlers = []  # clear handler list

    root_logger.addHandler(handler)
    root_logger.setLevel(logging.INFO)


def summarize_structure(structure):
    """Returns a dictionary with properties of the structure."""

    n_models = 0
    # per model
    n_chains = 0
    n_resids = 0
    n_ptmuts = 0  # DisorderedResidues
    n_icodes = 0
    n_aatoms = 0  # all atoms
    n_hatoms = 0  # hetatm
    n_altloc = 0

    for model in structure:
        n_models += 1

    for chain in model:  # all models should be the same anyway
        n_chains += 1
        for resid in chain:
            is_ptmut = resid.is_disordered() == 2
            n_ptmuts += int(is_ptmut)

            if is_ptmut:
                children = resid.disordered_get_list()
            else:
                children = [resid]

            for child in children:  # silly, but practical
                n_icodes += int(child.id[2] != ' ')
                n_resids += 1

                for atom in child.get_unpacked_list():
                    n_hatoms += int(atom.parent.id[0] != ' ')
                    n_aatoms += 1
                    n_altloc += int(bool(atom.altloc.strip()))

    return {
        'models': n_models,
        'chains': n_chains,
        'residues': n_resids,
        'res-icode': n_icodes,
        'res-ptmut': n_ptmuts,
        'all-atoms': n_aatoms,
        'het-atoms': n_hatoms,
        'altlocs': n_altloc
    }


def test_element_assignment(structure):
    for residue in structure.get_residues():
        for atom in residue.get_unpacked_list():
            # Test element assignment
            og_elem = atom.element.strip()
            atom._assign_element(element=None)  # force reassignment
            assert og_element and og_element == atom.element, \
                f'Element mismatch: "{og_elem}" != "{atom.element}"'


def main():

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        'infmt',
        choices=['pdb', 'mmcif'],
        help='File format of input files.'
    )
    ap.add_argument(
        'folder',
        type=pathlib.Path,
        help='Top-level folder with input files'
    )
    ap.add_argument(
        '--no-continue',
        action='store_true',
        default=False,
        help='Parses all input files, ignoring existing results.'
    )
    ap.add_argument(
        '--strict',
        action='store_true',
        default=False,
        help='Parse with PDBParser PERMISSIVE=0'
    )
    args = ap.parse_args()

    # Setup logging
    setup_logging()

    permissive_bool = not args.strict

    if args.infmt == 'pdb':
        parser = PDBParser(PERMISSIVE=permissive_bool, QUIET=1)
        writer = PDBIO()
    elif args.infmt == 'mmcif':
        parser = MMCIFParser(QUIET=1)
        writer = MMCIFIO()

    flist = sorted(args.folder.rglob('*.gz'))

    xmllist = sorted(args.folder.rglob('*.xml'))
    if not args.no_continue and xmllist:
        logging.info(f'Found {len(xmllist)} existing result files')
        xmlset = {
            f.stem: f for f in xmllist
        }
        fset = {f.stem: f for f in flist}
        remainder = set(fset.keys()) - set(xmlset.keys())
        logging.info(f'Resuming benchmark: {len(remainder)} files left')
        flist = sorted(fset[f] for f in remainder)
    else:
        logging.info(f'Found {len(flist)} files')

    n_digits = len(str(len(flist)))  # for fmting

    for idx, fpath in enumerate(flist, start=1):
        try:
            # Parse
            with gzip.open(fpath, mode='rt') as handle:
                t0 = time.time()
                s = parser.get_structure(fpath.name, handle)
                t1 = time.time()

                read_time = t1 - t0

                data = summarize_structure(s)

            # Write
            writer.set_structure(s)
            t0 = time.time()
            writer.save('io.temp')
            t1 = time.time()
            write_time = t1 - t0

            # Round-trip
            s2 = parser.get_structure('new', 'io.temp')
            data2 = summarize_structure(s2)

            assert data == data2, f'Summaries differ: {data} != {data2}'

        except Exception as err:
            with fpath.with_suffix('.failed').open('w') as f:
                print(err, file=f)
                print(traceback.format_exc(), file=f)

            status = 'failed'

        else:

            # Write XML file with numbers
            root = Element('structure')
            root.set('path', fpath.name)
            root.set('parse_time', f'{read_time:5.3f}')
            root.set('write_time', f'{write_time:5.3f}')

            for key, value in data.items():
                child = SubElement(root, key)
                child.text = str(value)

            # Reparse for pretty print
            xml = minidom.parseString(tostring(root, 'utf-8'))

            # Write to file
            with fpath.with_suffix('.xml').open('w') as f:
                f.write(xml.toprettyxml(indent='  '))

            # Clear XML memory
            root.clear()
            xml.unlink()
            del root, xml

            status = 'ok'

        finally:
            try:
                os.remove('io.temp')
            except Exception:
                pass

            memusage = psutil.virtual_memory().percent

            logging.info(
                f'{idx:>{n_digits}d}/{len(flist)} {fpath.parent.name}/{fpath.name}: {status} | mem% = {memusage}',
            )  # to check for leaks

if __name__ == '__main__':
    main()
