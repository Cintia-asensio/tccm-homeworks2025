#hw\_act

Program to compute closed-shell Hartree–Fock (HF) and MP2 energies from TREXIO (.h5) files.

## What it does
Given a TREXIO file containing molecular orbital integrals in the MO basis, the program reads the required data and computes the closed-shell Hartree–Fock total energy and the MP2 correlation and total energies.

The program reads the following data from the TREXIO file:
- nuclear repulsion energy (E\_nn)
- number of molecular orbitals (mo\_num)
- number of occupied orbitals (nocc, closed-shell)
- one-electron integrals (h\_pq)
- two-electron integrals (pq|rs)
- molecular orbital energies (eps\_p)

It then computes:
- Hartree–Fock total energy (E\_HF\_total)
- MP2 correlation energy (E\_MP2\_corr)
- MP2 total energy (E\_MP2\_total)

## Directory structure
src/        source code (main.c)
tests/      simple tests and reference checks
data/       input TREXIO .h5 files (not tracked by git)
results/    output results table (results.csv)
run\_all.sh  helper script to run all .h5 files
Makefile    build and test targets
INSTALL.md  installation and usage instructions
LICENSE     license information
AUTHORS     authorship information
README.md   this file

## Quick start
To build the program, run:
make clean
make

To run the program on a TREXIO file:
./main path/to/file.h5

Example:
./main data/h2o.h5

## Example output
File:   data/h2o.h5
E\_nn:   9.1949655588
mo\_num: 24
nocc:   5
E\_HF\_elec:  -85.2217642670
E\_HF\_total: -76.0267987082
E\_MP2\_corr: -0.2039599741
E\_MP2\_total: -76.2307586823

## Notes
- MP2 is computed using the closed-shell canonical expression with spatial orbitals.
- Two-electron integrals are read in sparse format and expanded using 8-fold symmetry.
- The current implementation builds a dense ERI tensor (mo\_num^4), which is suitable for the small test systems used in this project.

## Requirements
- C compiler supporting C99
- TREXIO library (headers and shared library)
- TREXIO .h5 file containing molecular orbital integrals

## License
See the LICENSE file.

## Authors
See the AUTHORS file.

