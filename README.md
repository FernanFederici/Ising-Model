This repository contains a python script which performs monte carlo simulations of a 2-dimensional ising model (assumed square lattice). The lattice energy, spin pair correlation, and the average spin of the lattice. As the simulation progresses, the user will be prompted with the efficiency of the monte carlo moves.
 
 This script allows the user to specify the side length of the lattice, the temperature of the lattice, the coupling constant between the spins, and the number of monte carlo moves to be implemented. Currently these values need to be set in the script, but it could be easily updated to be passed at the command line using argparse functionality. Something to do at another time!

## Potts Model C++ Code
This repository also contains a three-state Potts Model implemented with C++, which allows for an external field to be set along any of the three dimensions.

You will need a C++ compiler in order to turn potts.cpp and potts.hpp into an executable, most commonly available on Macs is the g++ compiler.

g++ potts.cpp -o pottsExecutable

./pottsExecutable

When you run the simulation you will be promted to input parameters for the simulation from the command line. These include, the temperature, coupling constants, external field diretion and magnitude, as well as the number of spins on the side of the lattice (assuming a square lattice), and the number of Monte Carlo sweeps of the lattice to perform.

This simulation generates a couple of files when executed, including potts.out and latticeConfig.out. potts.out contains the total energy of the lattice as a function of number of Monte Carlo sweeps, as well as the three components of magnetization. latticeConfig.out contains the resulting configuration of the lattice at the end of the simulation.
