#!/usr/bin/env python

"""
Script to benchmark Bio.PDB PDBParser on a collection of PDB files.

Parses all PDB files in a folder (expected gziped by the way),
writes them out, and parses them again.

Outputs XML files with some numbers on the parsed Structure
objects for later comparison with other datasets/references.
"""

import argparse
import gzip
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

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO


def main():

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        'folder',
        type=pathlib.Path,
        help='Top-level folder with PDB files'
    )
    ap.add_argument(
        '--strict',
        action='store_true',
        default=False,
        help='Parse with PERMISSIVE=0'
    )
    args = ap.parse_args()

    permissive_bool = not args.strict
    parser = PDBParser(PERMISSIVE=permissive_bool, QUIET=1)
    writer = PDBIO()

    pdblist = sorted(args.folder.rglob('*.ent.gz'))
    print(f'Found {len(pdblist)} files')

    for idx, pdbf in enumerate(pdblist, start=1):
        print(f'Parsing file: {pdbf.parent.name}/{pdbf.name}')
        try:
            # Parse
            with gzip.open(pdbf, mode='rt') as handle:
                mem0 = psutil.virtual_memory().percent
                t0 = time.time()
                s = parser.get_structure(pdbf.name, handle)
                t1 = time.time()
                mem1 = psutil.virtual_memory().percent

                read_time = t1 - t0
                read_memu = mem1 - mem0

            # Write
            writer.set_structure(s)
            t0 = time.time()
            writer.save('temp.pdb')
            t1 = time.time()
            write_time = t1 - t0

            # Round-trip
            s2 = parser.get_structure('new', 'temp.pdb')

            assert len(list(s.get_models())) == len(list(s2.get_models()))
            assert len(list(s.get_chains())) == len(list(s2.get_chains()))
            assert len(list(s.get_residues())) == len(list(s2.get_residues()))
            assert len(list(s.get_atoms())) == len(list(s2.get_atoms()))

        except Exception as err:
            with pdbf.with_suffix('.failed').open('w') as f:
                print(err, file=f)
                print(traceback.format_exc(), file=f)
        else:
            n_models = len(list(s.get_models()))
            n_chains = len(list(s.get_chains()))
            n_resids = len(list(s.get_residues()))
            n_atoms = len(list(s.get_atoms()))

            # Write XML file with numbers
            root = Element('structure')
            root.set('path', pdbf.name)
            root.set('parse_time', f'{read_time:5.3f}')
            root.set('write_time', f'{write_time:5.3f}')
            root.set('memory_usage', f'{read_memu}')

            for level in ('models', 'chains', 'residues', 'atoms'):
                child = SubElement(root, level)
                prop = getattr(s, f'get_{level}')
                child.text = str(len(list(prop())))

            # Reparse for pretty print
            xml = minidom.parseString(tostring(root, 'utf-8'))

            # Write to file
            with pdbf.with_suffix('.results.xml').open('w') as f:
                f.write(xml.toprettyxml(indent='  '))

            # Clear XML memory
            root.clear()
            xml.unlink()
            del root, xml

        finally:
            os.remove('temp.pdb')
            memusage = psutil.virtual_memory().percent
            print(
                f'Read {idx} out of {len(pdblist)} files. Memory %: {memusage}',
                flush=True
            )  # to check for leaks

if __name__ == '__main__':
    main()