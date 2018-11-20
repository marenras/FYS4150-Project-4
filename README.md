# FYS4150 - Computational physics: Project 4

In this respotory you will find codes used for project 4 in FYS4150.

## Description of programs

**Not MPI**

Compiled by writing `make` in command line. Can be runned in two different ways. If you want the program to run for several different temperature, the following variables has to be provided in the command line:  
* Lattice dimension
* Number of MC cycles 
* Initial temperature 
* Final temperature 
* Temperature step 
* Filename for output data

If you want the program to only run for one temperature, the following variables has to be provided in the command line:
* Lattice dimension 
* Number of MC cycles 
* Temperature 
* Filename for output data  


**MPI**

Compiled be writing `mpic++ -o run.x main_MPI.cpp ising.cpp` in command line. Runned by writing `mpirun ./run.x` and the following variables in the command line:
* Lattice dimension
* Number of MC cycles 
* Initial temperature 
* Final temperature 
* Temperature step 
* Filename for output data

It is also possible to add for example `-n 3` after `mpirun` if you want the program to run on three cores. If this is not provided, it will run on the number of cores given from the default host file.
