This repository contains the code for 2 programs, both will be created by running make.

The first program called isingnopar.exe is a simple program for running the metropolis algorithm without paralellization, primarily for the porpouse of debuging
the agrotim independent of MPI. 

The second program ising.exe is a parallellized version of isingnopar.exe.

Both programs take 3 inputs fisrtly the filename of the output file, secondly the lattice size and finaly the number of mc cycles. 
Note that the lattice will allways be square and inputing the integer n will relsult in an nxn lattice. For the mc cycles inputing the integer n results in 10^n cycles.

ising.exe is currently setup for 
