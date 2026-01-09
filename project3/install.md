How to install the dynamics program?

This program has been written in fortran90 and these is a compile language so in order to use the program we have to previous compile it.

The program is called molecular_dynamics.f90 

To compile the program you must write in the terminal gfortran molecular_dynamics.f90 -o molecular_dynamics 

Your input file must be called inp.txt and also must be in the same folder that the molecular_dynamics.f90 file

The format of this input file must be the first line containing the number of atoms and the subsequent line contains the x, y and z coordinates in Anstrong units followed by the mass of an atom.

In order to run the program you pust write on the terminal in the same folder ./molecular_dynamics

When the program finished calculating you will obtain a traj.xyz file that is the trajectory file.
