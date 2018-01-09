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

  ofstream energyFile;
  energyFile.open("potts.out");
  energyFile << "#iCycle \t Energy \t Magnetization[0] \t Magnetization[1] \t Magnetization[2] \n";
  
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
  int spinPreFlip = 0;
  int spinPostFlip = 0;
  
  float deltaE = 0.0;
  float deltaS = 0.0;

  // oldMagnetization contains the net # of spins aligned with each direction
  myLattice.calcMagnetization(oldMagnetization);
  
  energyFile << "0\t" + std::to_string(oldEnergy) + "\t" + 
    std::to_string(oldMagnetization[0]/nSpins) + "\t" +
    std::to_string(oldMagnetization[1]/nSpins) + "\t" +
    std::to_string(oldMagnetization[2]/nSpins) + "\n";
  
  std::cout << "Energy of lattice at the beginning = " << oldEnergy << endl;
  std::cout << "Magnetization of lattice at the beginning = "
	    << oldMagnetization[0]/nSpins << ","
	    << oldMagnetization[1]/nSpins << ","
	    << oldMagnetization[2]/nSpins << endl;
  
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
	spinPreFlip = myLattice.getSpin(idex,jdex);
	
	myLattice.flipSpin(idex,jdex);
	
	postFlipEnergy = myLattice.calcDifferenceInEnergy(idex,jdex);
	spinPostFlip = myLattice.getSpin(idex,jdex);
	
	deltaE = postFlipEnergy - preFlipEnergy;
	newEnergy = oldEnergy + deltaE;

	  
	if (newEnergy > oldEnergy){
	  double randNum = (double)rand()/(double)RAND_MAX;
	  // compare new energy to a thermal distribution
	  if (randNum <= exp(-1.0/myLattice.getTemperature())*(newEnergy - oldEnergy)) {
	    acceptedMoves = acceptedMoves + 1;
	    oldEnergy = newEnergy;
	    oldMagnetization[spinPreFlip - 1] -= 1;
	    oldMagnetization[spinPostFlip - 1] += 1;

	  }
	  // new energy is too large when compared to a thermal distribution
	  // flip spin again to go back to old configuration.
	  else {
	    myLattice.revertSpin(idex,jdex, spinPreFlip);
	  }
	}
	
	// if the new energy configuration is lower than old, keep the move
	else {
	  acceptedMoves = acceptedMoves + 1;
	  oldEnergy = newEnergy;
	  oldMagnetization[spinPreFlip - 1] -= 1;
	  oldMagnetization[spinPostFlip - 1] += 1;	    

	}
	
      } // idex loop
    } // jdex loop

    energyFile << std::to_string(iCycle+1) + "\t" +
      std::to_string(oldEnergy) + "\t" + 
      std::to_string(oldMagnetization[0]/nSpins) + "\t" + 
      std::to_string(oldMagnetization[1]/nSpins) + "\t" + 
      std::to_string(oldMagnetization[2]/nSpins) + "\n";


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


   ofstream latticeOutFile;
   latticeOutFile.open("latticeConfig.out");
   latticeOutFile << "#Configuration of the lattice after simulation \n";
   for (int idex = 0; idex < myLattice.getSideLength(); idex++){
     latticeOutFile <<  "   ";
     for (int jdex = 0; jdex < myLattice.getSideLength(); jdex++){	
       latticeOutFile << std::to_string(myLattice.getSpin(idex,jdex)) + "  ";
     }
     latticeOutFile << "\n";
   }
  
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
