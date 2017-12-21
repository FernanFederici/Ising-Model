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
  float externFieldDirect;
  
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

  int getExternFieldDirect(void);
  void setExternFieldDirect(float newFieldDirect);
  
  // Utility functions of the lattice.
  void initializeLattice(void);
  void printLattice(void);
  float calcTotalEnergy(void);
  void calcMagnetization(float magnetization[3]);
  int getRandomCoord(void);
  float flipSpin(int idex, int jdex);
  void revertSpin(int idex, int jdex, float savedState);
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

int Lattice::getExternFieldDirect(void){
  return externFieldDirect;
}

void Lattice::setExternFieldDirect(float newFieldDirect){
  externFieldDirect = newFieldDirect;
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

  string localExternalFieldDirect;
  std::cout << "Please enter the external field direction (1.0, 2.0, 3.0):" << endl;
  std::cin >> localExternalFieldDirect;
  float flocalExternalFieldDirect = std::stof(localExternalFieldDirect);
  setExternFieldDirect(flocalExternalFieldDirect);

  string localExternalField;
  std::cout << "Please enter the external field:" << endl;
  std::cin >> localExternalField;
  float flocalExternalField = std::stof(localExternalField);
  setExternField(flocalExternalField);
  
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

    H = -J * Sum( kron[si,sj] ) - h * Sum( kron[hi,si] )

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
	if (matrix[i][j] == matrix[idex][j]) 
	  energySum = energySum - 1.0;
      }
      else if ( i-1 >= 0) {
	idex = i - 1;
	if (matrix[i][j] == matrix[idex][j]) 
	  energySum = energySum - 1.0;
      }
      
      // (i,j) <-> (i+1,j)
      if ( i+1 == sideLength) {
	idex = i + 1 - sideLength;
	if (matrix[i][j] == matrix[idex][j]) 
	  energySum = energySum - 1.0;
      }
      else if (i+1 < sideLength) {
	idex = i + 1;
	if (matrix[i][j] == matrix[i][j]) 
	  energySum = energySum - 1.0;
      }
      
      // (i,j) <-> (i,j-1)
      if (j-1 < 0) {
	jdex = j - 1 + sideLength;
	if (matrix[i][j] == matrix[i][jdex])
	  energySum = energySum - 1.0;
      }
      else if (j-1 >= 0) {
	jdex = j - 1;
	if (matrix[i][j] == matrix[i][jdex])
	  energySum = energySum - 1.0;
      }
      
      // (i,j) <-> (i,j+1)
      if (j+1 == sideLength) {
	jdex = j + 1 - sideLength;
	if (matrix[i][j] == matrix[i][jdex])
	  energySum = energySum - 1.0;
      }
      else if (j+1 < sideLength){
	jdex = j;
	if (matrix[i][j] == matrix[i][jdex]) 
	  energySum = energySum - 1.0;
      }

      if (matrix[i][j] == externFieldDirect)
	energySum = energySum - externField;
      
    }
  }
  energySum = energySum * coupling;
  return energySum;
}

void Lattice::calcMagnetization(float magnetization[3]) {
  float netA = 0.0;
  float netB = 0.0;
  float netC = 0.0;
  
  for (int i = 0; i < sideLength; i++) {
    for (int j = 0; j < sideLength; j++) {
      
      if (matrix[i][j] == 1.0)
	netA = netA + 1.0;
      else if (matrix[i][j] == 2.0)
	netB = netB + 1.0;
      else if (matrix[i][j] == 3.0)
	netC = netC + 1.0;
      
    }
  }
  
  netA = netA / pow(sideLength, 2.0);
  netB = netB / pow(sideLength, 2.0);
  netC = netC / pow(sideLength, 2.0);

  magnetization[0] = netA;
  magnetization[1] = netB;
  magnetization[2] = netC;

}

int Lattice::getRandomCoord(void) {
  double val = (double)rand()/(double)RAND_MAX;
  return (val*sideLength);
}

float Lattice::flipSpin(int idex, int jdex) {
  double val = (double)rand()/(double)RAND_MAX;
  float savedState = matrix[idex][jdex];
  
  if (matrix[idex][jdex] == 1.0) {
    if (val < 0.5) {
      matrix[idex][jdex] = 2.0;
      return savedState;
    }
    else {
      matrix[idex][jdex] = 3.0;
      return savedState;
    }
  }
    
  if (matrix[idex][jdex] == 2.0) {
    if (val < 0.5) {
      matrix[idex][jdex] = 3.0;
      return savedState;
    }
    else {
      matrix[idex][jdex] = 1.0;
      return savedState;
    }
  }

  if (matrix[idex][jdex] == 3.0) {
    if (val < 0.5) {
      matrix[idex][jdex] = 1.0;
      return savedState;
    }
    else {
      matrix[idex][jdex] = 2.0;
      return savedState;
    }
  }
  return savedState;
}

void Lattice::revertSpin(int idex, int jdex, float savedState) {
  matrix[idex][jdex] = savedState;
}

float Lattice::calcDifferenceInEnergy(int idex, int jdex) {
  float randSpinE = 0.0;
  int i = 0;
  int j = 0;
  /*
    Here we will only consider neighbors of the flipped spin for the i,j pairs

    H = -J * Sum( kron[si,sj] ) - h * Sum( kron[hi, si] )

    There are four pairs of energies that need to be computed for each i,j pair.
    for a given i,j, there are 
    (i,j) <-> (i-1,j)     (i,j) <-> (i,j-1)
    (i,j) <-> (i+1,j)     (i,j) <-> (i,j+1)

  */

   // (i,j) <-> (i-1,j)
  if ( idex-1 < 0){
    i = idex - 1 + sideLength;
    if (matrix[idex][jdex] == matrix[i][jdex])
      randSpinE = randSpinE - 1.0;
  }
  else if (idex-1 >= 0) {
    if (matrix[idex][jdex] == matrix[idex-1][jdex])
      randSpinE = randSpinE - 1.0;
  }

  // (i,j) <-> (i+1,j)
  if (idex+1 == sideLength) {
    i = idex + 1 - sideLength;
    if (matrix[idex][jdex] == matrix[i][jdex])
      randSpinE = randSpinE - 1.0;
  }
  else if (idex+1 < sideLength) {
    if (matrix[idex][jdex] == matrix[idex+1][jdex])
      randSpinE = randSpinE - 1.0;
  }

  // (i,j) <-> (i,j-1)
  if ( jdex-1 < 0){
    j = jdex - 1 + sideLength;
    if (matrix[idex][jdex] == matrix[idex][j])
      randSpinE = randSpinE - 1.0;
  }
  else if (jdex-1 >= 0) {
    if (matrix[idex][jdex] == matrix[idex][jdex-1])
      randSpinE = randSpinE - 1.0;
  }

  // (i,j) <-> (i,j+1)
  if (jdex+1 == sideLength) {
    j = jdex + 1 - sideLength;
    if (matrix[idex][jdex] == matrix[idex][j])
      randSpinE = randSpinE - 1.0;
  }
  else if (jdex+1 < sideLength) {
    if (matrix[idex][jdex] == matrix[idex][jdex+1])
      randSpinE = randSpinE - 1.0;
  }


  
  if (matrix[idex][jdex] == externFieldDirect)
    randSpinE = randSpinE - externField;
  
  
  
  randSpinE = randSpinE * 2.0 * coupling;
  //std::cout << "randSpinE = " << randSpinE << endl;
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


