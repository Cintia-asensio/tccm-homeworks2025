# Installation and usage

This document describes how to install the required dependencies, build the program, and run it on a TREXIO file.

## Requirements
- C compiler supporting C99 (gcc or clang)
- TREXIO library
- HDF5 library (required by TREXIO)
- make

## Installing TREXIO
TREXIO can be installed from source or using a package manager, depending on the system. The instructions below describe a simple installation from source.
```
git clone https://github.com/TREX-CoE/trexio.git
cd trexio
cmake -B build
cmake --build build
cmake --install build
```
## Building
To build the program from the repository root:
```
make clean
make
```
After building, an executable named main should be available.

## Running the program
To run the program, provide a TREXIO .h5 file as command-line argument.
```
./main path/to/file.h5
```
Example:
```
./main data/h2o.h5
```
