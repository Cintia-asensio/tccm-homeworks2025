## A molecular Dynamics Code

### What it does
Molecular dynamics simulates the movement of atoms based on their initial positions and velocities. This program read force field parameters and the initial positions of atoms from a input file called "inp.txt".

The program reads the following data form the inp.txt file:
- For the first line it reads the number of atoms
- In the following lines it reads the x, y, z coordinates in Anstrong followed by the mass of an atom.

Then the program computes:
- The internuclear distance between each pair of atoms using euclidean formula.
- The Leannard-Jones potential
- The total kinetic energy of the system and the total energy of the system.
- The aceleration vector for each atom and stores it in a double precission array.
- Finnaly it implements the molecular dynamics

## Directory structure:
In the src directory there is the source code called "molecular_dynamics.f90" and a "inp.txt" file to have a example. In this way the user can see that it is mandatory to have this input file in the same directory that the code.

In the test directory you can find some test to ensure that the program works at it is expected.
 
## Quick start 
In order to build the program write
gfortran -o dynamics molecular_dynamics.f90

To run the program execute
./dynamics

## Output
The output contains a line with the number of atoms, a comment line with the step number and for each atom a line containing the atomic symbol followed by its x, y and z coordinates in Anstrong. This is written for every 10 steps of the dynamics.

## Notes
To compute this code there are some data values that if you want to change you must modify the molecular_dynamics.f90 file. This data values are:
- mass = 39.948 g/mol
- epsilon = 0.997 kJ/mol
- sigma = 3.405 Anstrongs
- time step = 0.02 picoseconds
- total number of steps = 1000

Atom coordinates will be read and written in Anstrong in input and output files, abd the code converts it into nanometers to have consistent units. The time is expresed in picoseconds.

## Requirements
- gfortran compiler

## License 
See the license file

## Authors 
See the authors file

