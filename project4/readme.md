This repository contains the code for 2 programs, both will be created by running make.

The first program called isingnopar.exe is a simple program for running the metropolis algorithm without paralellization, primarily for the porpouse of debuging
the agrotim independent of MPI. 

The second program ising.exe is a parallellized version of isingnopar.exe.

Both programs take 3 inputs firstly the filename of the output file, secondly the lattice size and finaly the number of mc cycles. 
Note that the lattice will allways be square and inputing the integer n will relsult in an nxn lattice. For the mc cycles inputing the integer n results in 10^n cycles.

ising.exe is currently setup for caclulating the energy, heatcapacity, magnetic moment and magnetic suseptibilty in the temprature range 2.0-2.3. in order to change the intilal/final temprature or the tempraure step these parametres must be changed in the code itself given in main.cpp.

the code for isingnopar.exe is given in nopar.cpp, it contains code for randomly initializing the lattice, looping over the number of mccycles and returing the number of accepted configurations or the energy and magnetic moment as a function of tim/mccycles. These features are currently commented out.
