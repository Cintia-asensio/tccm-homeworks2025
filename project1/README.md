# hw_act_mp2

Program to compute closed-shell Hartree–Fock (HF) and MP2 energies from TREXIO (`.h5`) files.

## What it does
Given a TREXIO file containing molecular orbital (MO) integrals in the MO basis, the program reads the required data and computes the closed-shell Hartree–Fock total energy and the MP2 correlation and total energies.

The program reads the following data from the TREXIO file:
- nuclear repulsion energy (`E_nn`)
- number of molecular orbitals (`mo_num`)
- number of occupied orbitals (`nocc`, closed-shell)
- one-electron integrals (`h_pq`)
- two-electron integrals (`<pq|rs>`)
- molecular orbital energies (`eps_p`)

It then computes:
- Hartree–Fock total energy (`E_HF_total`)
- MP2 correlation energy (`E_MP2_corr`)
- MP2 total energy (`E_MP2_total`)

## Directory structure
```
src/         source code (main.c)
tests/       tests and reference checks
data/        input TREXIO .h5 files (not tracked by git)
results/     output results table (results.csv)
run_all.sh   helper script to run all .h5 files
Makefile     build and test targets
INSTALL.md   installation and usage instructions
LICENSE      license information
AUTHORS      authorship information
README.md    this file
```

## Quick start

### Build
Compile the program using:
```bash
make clean
make
```

This produces the executable:
```
mp2_energy
```

### Run on a single file
Run the program on a single TREXIO file:
```bash
./mp2_energy path/to/file.h5
```

Example:
```bash
./mp2_energy data/h2o.h5
```

### Run tests (numerical validation)
To check that the implementation is correct, run:
```bash
make test
```

This command:
- builds `mp2_energy` if needed
- runs reference tests (currently on H2O)
- checks HF and MP2 energies against known reference values within a tolerance

### Run all molecules and generate results
To run the program on all `.h5` files in `data/` and generate a CSV table:
```bash
chmod +x run_all.sh
./run_all.sh
```

The results are written to:
```
results/results.csv
```

## Example output
```
File:   data/h2o.h5
E_nn:   9.1949655588
mo_num: 24
nocc:   5
E_HF_elec:  -85.2217642670
E_HF_total: -76.0267987082
E_MP2_corr: -0.2039599741
E_MP2_total: -76.2307586823
```

## Notes
- MP2 is computed using the closed-shell canonical expression with spatial orbitals.
- Two-electron integrals are read in sparse format and expanded using 8-fold symmetry.
- The current implementation builds a dense ERI tensor (`mo_num^4`), which is suitable for the small test systems used in this project.

## Requirements
- C compiler supporting C99
- TREXIO library (headers and shared library)
- TREXIO `.h5` file containing molecular orbital integrals

## License
See the LICENSE file.

## Authors
See the AUTHORS file.

