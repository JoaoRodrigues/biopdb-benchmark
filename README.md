# Benchmarking Bio.PDB

These scripts allow you to download and benchmark the PDB and mmCIF
parser in Biopython (Bio.PDB) on the entire PDB database.

The benchmark scripts will produce a *.result.xml file for each data
file with timings and memory usage, as well as brief details of the
chemical composition of the structure (number of chains, atoms, etc).

The parsers are validated by cross-comparison and independent comparison
to similar XML output files derived from PDBML files (see download_pdb.sh).

## Instructions
```bash

pip install psutil  # to track memory usage
python -m pip install git+https://github.com/biopython/biopython.git  # get latest revision
./download_pdb.sh  # takes a few hours

python benchmark_PDBParser.py pdb/ &> benchmark.pdb.out
python benchmark_MMCIFParser.py pdb/ &> benchmark.cif.out
python process_xmldb.py XML &> xmlprocess.out
```

## Questions?
Feel free to send me an email or open an issue here and tag me.
