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

  ofstream outputFile;
  outputFile.open("potts.out");
  outputFile << "#iCycle \t Energy \n";
  
  int nSpins = pow(myLattice.getSideLength(),2);
  int sampleTime = 1;
  int attemptedMoves = 0;
  int acceptedMoves = 0; 
  float oldEnergy = myLattice.calcTotalEnergy();
  float oldMagnetization [3] = {0.0, 0.0, 0.0};
  float newEnergy = 0.0;
  float newMagnetization = 0.0;
  float preFlipEnergy = 0.0;
  float postFlipEnergy = 0.0;
  float savedState = 0.0;

  float deltaE = 0.0;
  float deltaS = 0.0;

  myLattice.calcMagnetization(oldMagnetization);
   outputFile << "0\t" + std::to_string(oldEnergy) + "\n";
  std::cout << "Energy of lattice at the beginning = " << oldEnergy << endl;
  std::cout << "Magnetization of lattice at the beginning = "
	    << oldMagnetization[0] << "," << oldMagnetization[1] << ","
	    << oldMagnetization[2] << endl;
  
  myLattice.printLattice();
  
  //Loop over the designated number of MCSweeps
  for (int iCycle = 0; iCycle < nCycles; iCycle++) {
    
    
    //For each sweep, loop through entire lattice
    for (int idex = 0; idex < myLattice.getSideLength(); idex++){
      for (int jdex = 0; jdex < myLattice.getSideLength(); jdex++){	
	
	attemptedMoves = attemptedMoves + 1;

	/*
	  With the Potts model and # states > 2, we must calculate the deltaE
	  of the spin flip explicitly, due to the possible interactions.
	*/
	
	preFlipEnergy = myLattice.calcDifferenceInEnergy(idex,jdex);
	savedState = myLattice.flipSpin(idex,jdex);
	postFlipEnergy = myLattice.calcDifferenceInEnergy(idex,jdex);
	
	deltaE = postFlipEnergy - preFlipEnergy;
	newEnergy = oldEnergy + deltaE;

	  
	if (newEnergy > oldEnergy){
	  double randNum = (double)rand()/(double)RAND_MAX;
	  // compare new energy to a thermal distribution
	  if (randNum <= exp(-1.0/myLattice.getTemperature())*(newEnergy - oldEnergy)) {
	    acceptedMoves = acceptedMoves + 1;
	    oldEnergy = newEnergy;

	  }
	  // new energy is too large when compared to a thermal distribution
	  // flip spin again to go back to old configuration.
	  else {
	    myLattice.revertSpin(idex,jdex, savedState);
	  }
	}
	
	// if the new energy configuration is lower than old, keep the move
	else {
	  acceptedMoves = acceptedMoves + 1;
	  oldEnergy = newEnergy;

	}
	
      } // idex loop
    } // jdex loop

    outputFile << std::to_string(iCycle+1) + "\t" + std::to_string(oldEnergy) + "\n";
  }// iCycle loop
  myLattice.calcMagnetization(oldMagnetization);
  std::cout << "Magnetization of lattice at the end = "
	    << oldMagnetization[0] << "," << oldMagnetization[1] << ","
	    << oldMagnetization[2] << endl;
  
  std::cout << "Printing out lattice after all MC steps are completed" << endl;
  myLattice.printLattice();

  std::cout << "Energy at end = " << oldEnergy << endl;

  std::cout << "attemptedMoves = " << attemptedMoves << endl;
  std::cout << "acceptedMoves = " << acceptedMoves << endl;
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
