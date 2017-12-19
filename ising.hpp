#include <vector>
#include <string>
#include <iostream>
#include <math.h>
using namespace std;

class Lattice {
public:
  float temperature;
  float coupling;
  int sideLength;
  int numSpins;
  float spinMoment;
  float externField;
  
  vector< vector<float> > matrix;
  vector<float> aveSpin;
  vector<float> latticeEnergy;
  vector<float> spinPairCorr;

  // Member functions declared.
  int getTemperature(void);
  void setTemperature(float newTemp);

  float getCoupling(void);
  void setCoupling(float newCoupling);
  
  int getSideLength(void);
  void setSideLength(int newSideLength);

  int getExternField(void);
  void setExternField(float newField);
  
  // Utility functions of the lattice.
  void initializeLattice(void);
  void printLattice(void);
  float calcTotalEnergy(void);
  float calcMagnetization(void);
  int getRandomCoord(void);
  void flipSpin(int idex, int jdex);
  float calcDifferenceInEnergy(int idex, int jdex);
  float calcDifferenceInMagnetization(int idex, int jdex);
  
  void accumulateStats(float energy, float magnetization);

  
};




// Member function definitions
int Lattice::getTemperature(void) {
  return temperature;
}

void Lattice::setTemperature(float newTemp) {
  if (newTemp >= 0.0) 
    temperature = newTemp;
  
  else 
    std::cout << "Temperature must be >= 0.0" << endl;
}

float Lattice::getCoupling(void) {
  return coupling;
}

void Lattice::setCoupling(float newCoupling) {
  coupling = newCoupling;
}

int Lattice::getSideLength(void) {
  return sideLength;
}

void Lattice::setSideLength(int newSideLength) {
  if (newSideLength > 0) {
    sideLength = newSideLength;
  }
  else
    std::cout << "Side length must be > 0." << endl; 
}

int Lattice::getExternField(void){
  return externField;
}

void Lattice::setExternField(float newField){
  externField = newField;
}


// Utility Functions
/*
initializeLattice will ask the user to set the temperature, coupling 
constant, and side length (num of spins to a side).
The function then sets the parameters and builds the lattice.
*/
void Lattice::initializeLattice(void) {
  string localTemp;
  std::cout << "Please enter the temperature:" << endl;
  std::cin >> localTemp;
  float flocalTemp = std::stof(localTemp);
  setTemperature(flocalTemp);

  string localCoupling;
  std::cout << "Please enter the coupling constant:" << endl;
  std::cin >> localCoupling;
  float flocalCoupling = std::stof(localCoupling);
  setCoupling(flocalCoupling);
  
  string localSideLength;
  std::cout << "Please enter integer N spins for NxN lattice:" << endl;
  std::cin >> localSideLength;
  int ilocalSideLength = std::atoi((localSideLength).c_str());
  setSideLength(ilocalSideLength);

  // Next we will initialize the matrix (lattice) to all spin up
  for (int i = 0; i < ilocalSideLength; i++) {
    vector<float> row;
    for (int j = 0; j < ilocalSideLength; j++) {
      row.push_back(1.0);
    }
    matrix.push_back(row);
  }

  // Print out values that we have set
  //std::cout << "Temp set to " << getTemperature() << endl;
  //std::cout << "Coupling set to " << getCoupling() << endl;
  //std::cout << "SideLength set to " << getSideLength() << endl;
 
}


void Lattice::printLattice(void) {
  for (int i = 0; i < sideLength; i++) {
    for (int j = 0; j < sideLength; j++) {
      std::cout << matrix[i][j] << "\t";
    }
    std::cout << "\n";
  }
}


float Lattice::calcTotalEnergy(void) {
  float energySum = 0.0;
  int idex = 0;
  int jdex = 0;
  /*
    There are four pairs of energies that need to be computed for each i,j pair.
    for a given i,j, there are 
    (i,j) <-> (i-1,j)     (i,j) <-> (i,j-1)
    (i,j) <-> (i+1,j)     (i,j) <-> (i,j+1)
   */

  
  for (int i = 0; i < sideLength; i++) {
    for (int j = 0; j < sideLength; j++) {
      
      // (i,j) <-> (i-1,j)
      if ( i-1 < 0) {
	idex = i - 1 + sideLength;
	energySum = energySum - (matrix[i][j] * matrix[idex][j]);
      }
      else if ( i-1 >= 0) {
	  idex = i - 1;
	  energySum = energySum - (matrix[i][j] * matrix[idex][j]);
	}

      // (i,j) <-> (i+1,j)
      if ( i+1 == sideLength) {
	idex = i + 1 - sideLength;
	energySum = energySum - (matrix[i][j] * matrix[idex][j]);
      }
      else if (i+1 < sideLength) {
	  idex = i + 1;
	  energySum = energySum - (matrix[i][j] * matrix[i][j]);
	}

      // (i,j) <-> (i,j-1)
      if (j-1 < 0) {
	jdex = j - 1 + sideLength;
	energySum = energySum - (matrix[i][j] * matrix[i][jdex]);
      }
      else if (j-1 >= 0) {
	  jdex = j - 1;
	  energySum = energySum - (matrix[i][j] * matrix[i][jdex]);
	}

      // (i,j) <-> (i,j+1)
      if (j+1 == sideLength) {
	jdex = j + 1 - sideLength;
	energySum = energySum - (matrix[i][j] * matrix[i][jdex]);
      }
      else if (j+1 < sideLength){
	jdex = j;
	energySum = energySum - (matrix[i][j] * matrix[i][jdex]);
      }
	
    }
  }
  energySum = energySum * (coupling / 2.0);
  return energySum;
}

float Lattice::calcMagnetization(void) {
  float magnetization = 0.0;
  for (int i = 0; i < sideLength; i++) {
    for (int j = 0; j < sideLength; j++) {
      magnetization = magnetization + matrix[i][j];
    }
  }
  return (magnetization / pow(sideLength, 2.0));
}

int Lattice::getRandomCoord(void) {
  double val = (double)rand()/(double)RAND_MAX;
  return (val*sideLength);
}

void Lattice::flipSpin(int idex, int jdex) {
  matrix[idex][jdex] = -1.0 * matrix[idex][jdex];
  std::cout << matrix[idex][jdex] << endl;
}


float Lattice::calcDifferenceInEnergy(int idex, int jdex) {
  float randSpinE = 0.0;
  int i = 0;
  int j = 0;
  /*
    This time we will have to be more careful with neighboring interactions.
  */

  if ( idex-1 < 0){
    i = idex - 1 + sideLength;
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[i][jdex]);
  }
  else if (idex-1 >= 0) {
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex-1][jdex]);
  }

  if (idex+1 == sideLength) {
    i = idex + 1 - sideLength;
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[i][jdex]);
  }
  else if (idex+1 < sideLength) {
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex+1][jdex]);
  }

   if ( jdex-1 < 0){
    j = jdex - 1 + sideLength;
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex][j]);
  }
  else if (jdex-1 >= 0) {
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex][jdex-1]);
  }

  if (jdex+1 == sideLength) {
    j = jdex + 1 - sideLength;
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex][j]);
  }
  else if (jdex+1 < sideLength) {
    randSpinE = randSpinE - (matrix[idex][jdex] * matrix[idex][jdex+1]);
  }

  randSpinE = randSpinE * 2.0 * coupling;
  return randSpinE;
}


float Lattice::calcDifferenceInMagnetization(int idex, int jdex) {
  float newMagnetization;
  newMagnetization = 2.0 * matrix[idex][jdex] / pow(sideLength,2);

  return newMagnetization;
}
 

void Lattice::accumulateStats(float energy, float magnetization) {
  aveSpin.push_back(magnetization);
  latticeEnergy.push_back(energy); 
}


