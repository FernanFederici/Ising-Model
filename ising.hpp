#include <vector>
#include <string>
#include <iostream>
using namespace std;

class Lattice {
public:
  float temperature;
  float coupling;
  int sideLength;
  int numSpins;
  float spinMoment;

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

  void initializeLattice(void);

  
};




// Member function definitions
int Lattice::getTemperature(void) {
  return temperature;
}

void Lattice::setTemperature(float newTemp) {
  if (newTemp >= 0.0) {
    temperature = newTemp;
  }
  else {
    std::cout << "Temperature must be >= 0.0" << endl;
  }
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
  sideLength = newSideLength;
}



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

  
}



