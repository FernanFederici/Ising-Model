/*
  This program will simulate a two-dimensional potts-model.
 */

#include "ising.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "string"
using namespace std;


void readFile(string fileName) {
  ifstream myReadFile;
  myReadFile.open(fileName);
}



void MonteCarloSim(Lattice myLattice,int nCycles) {
  int nSpins = pow(myLattice.getSideLength(),2);
  int sampleTime = 1;
  int attemptedMoves = 0;
  int acceptedMoves = 0; 
  float oldEnergy = myLattice.calcTotalEnergy();
  float oldMagnetization = myLattice.calcMagnetization();
  float newEnergy = 0.0;
  float newMagnetization = 0.0;

  float deltaE = 0.0;
  float deltaS = 0.0;

  std::cout << "Energy of lattice at the beginning = " << oldEnergy << endl;
  myLattice.printLattice();
  
  //Loop over the designated number of MCSweeps
  for (int iCycle = 0; iCycle < nCycles; iCycle++) {
    
    myLattice.accumulateStats(oldEnergy,oldMagnetization);
    
    //For each sweep, loop through entire lattice
    for (int idex = 0; idex < myLattice.getSideLength(); idex++){
      for (int jdex = 0; jdex < myLattice.getSideLength(); jdex++){	
	
	attemptedMoves = attemptedMoves + 1;

	myLattice.flipSpin(idex,jdex);
	deltaE = myLattice.calcDifferenceInEnergy(idex,jdex);
	deltaS = myLattice.calcDifferenceInMagnetization(idex,jdex);
	
	newEnergy = oldEnergy + deltaE;
	newMagnetization = oldMagnetization + deltaS;
	  
	if (newEnergy > oldEnergy){
	  double randNum = (double)rand()/(double)RAND_MAX;
	  // compare new energy to a thermal distribution
	  if (randNum <= exp(-1.0/myLattice.getTemperature())*(newEnergy - oldEnergy)) {
	    acceptedMoves = acceptedMoves + 1;
	    oldEnergy = newEnergy;
	    oldMagnetization = newMagnetization;
	  }
	  // new energy is too large when compared to a thermal distribution
	  // flip spin again to go back to old configuration.
	  else {
	    myLattice.flipSpin(idex,jdex);
	  }
	}
	// if the new energy configuration is lower than old, keep the move
	else {
	  acceptedMoves = acceptedMoves + 1;
	  oldEnergy = newEnergy;
	  oldMagnetization = newMagnetization;
	}


	//if (iCycle % (nCycles/1000) == 0) {
	  //myLattice.printLattice();
	  //std::cout << "energy of lattice = " << oldEnergy << endl;
	  //std::cout << "magnetization of lattice = " << oldMagnetization << endl;
	  //float percent = float(iCycle) / float(nCycles);
	  //std::cout << "%complete = " << percent * 100.0  << endl;
	//}
	
      } // idex loop
    } // jdex loop

  }// iCycle loop
  std::cout << "Printing out lattice after all MC steps are completed" << endl;
  myLattice.printLattice();

  std::cout << "Energy at end = " << oldEnergy << endl;
}


int main() {

  // string fileName;
  // std::cout << "Please enter a file name to be read in." << endl;
  // std::cin >> fileName;
  // readFile(fileName);

  
  Lattice myLattice;
  myLattice.initializeLattice();

  int nCycles;
  std::cout << "Please enter the number of Monte Carlo Sweeps:" << endl;
  std::cin >> nCycles;
  
  MonteCarloSim(myLattice, nCycles);
  
  return 0;
}
