import random
import math
#import matplotlib.pyplot as plt
import numpy as np

class Lattice:
    def __init__(self,sideLength):
        self.temperature = 1.0
        self.coupling = 1.0
        self.sideLength = sideLength
        self.numSpins = sideLength * sideLength
        self.matrix = [[math.pow(1.0,j) for j in range(sideLength)]for i in range(sideLength)]
        self.aveSpin = []
        self.latticeEnergy = []
        self.moment = 1.0
        self.spinPairCorr = []

    def getSideLength(self):
        return self.sideLength
    def getMatrix(self):
        return self.matrix
    def setTemperature(self,newTemperature):
        self.temperature = newTemperature    
    def getTemperature(self):
        return self.temperature
    def setCoupling(self,newCoupling):
        self.coupling = newCoupling
    def getCoupling(self):
        return self.coupling

    def printMatrix(self):
        for i in range(0, self.sideLength):
	  for j in range(0, self.sideLength):
	    print str(self.matrix[i][j]) + "\t",
	  print ""

    def calcTotalEnergy(self):
        energySum = 0.0
        for i in range(0,self.getSideLength()):
            for j in range(0,self.getSideLength()):
                energySum = energySum - self.getMatrix()[i][j] * self.getMatrix()[i-1][j]
                energySum = energySum - self.getMatrix()[i][j] * self.getMatrix()[i-self.getSideLength()+1][j]
                energySum = energySum - self.getMatrix()[i][j] * self.getMatrix()[i][j-1]
                energySum = energySum - self.getMatrix()[i][j] * self.getMatrix()[i][j-self.getSideLength()+1]
        energySum = energySum * (self.getCoupling() / 2.0)
        return energySum

    def calcMagnetization(self):
        magnetization = 0.0
        for i in range(0, self.getSideLength()):
            for j in range(0, self.getSideLength()):
                magnetization = magnetization + self.getMatrix()[i][j]
        return (magnetization / math.pow(self.getSideLength(),2))
    

    def spinPair(self,distance):
        spinpair = 0.0
        for i in range(0, self.getSideLength()):  
            for j in range(0, self.getSideLength()):
                if int(i + distance) >= int(self.getSideLength()):             
                    ineighbor = (i + distance) - self.getSideLength()     
                else: ineighbor = (i + distance)
                spinpair = spinpair + self.getMatrix()[i][j] * self.getMatrix()[ineighbor][j]
                if int(j + distance) >= int(self.getSideLength()):          
                    jneighbor = (j + distance) - self.getSideLength()         
                else: jneighbor = (j + distance)                           
                spinpair = spinpair + self.getMatrix()[i][j] * self.getMatrix()[i][jneighbor]                                                        
        normSpinPair = spinpair / (2*math.pow(self.getSideLength(),2))         
        return normSpinPair                                       



    def getRandomSpin(self):
        randomX = int(random.uniform(0,1)*self.getSideLength())
        randomY = int(random.uniform(0,1)*self.getSideLength())
        randomSpin = [randomX, randomY]
        return randomSpin

    def flipSpin(self,randSpin):
        self.getMatrix()[randSpin[0]][randSpin[1]] = (-1.0) * self.getMatrix()[randSpin[0]][randSpin[1]]

    def calcDifferenceInEnergy(self,randSpin):
        randSpinE = 0.0
        randSpinE = randSpinE - self.getMatrix()[randSpin[0]][randSpin[1]] * self.getMatrix()[randSpin[0]-1][randSpin[1]]
        randSpinE = randSpinE - self.getMatrix()[randSpin[0]][randSpin[1]] * self.getMatrix()[randSpin[0]-self.getSideLength()+1][randSpin[1]]
        randSpinE = randSpinE - self.getMatrix()[randSpin[0]][randSpin[1]] * self.getMatrix()[randSpin[0]][randSpin[1]-1]
        randSpinE = randSpinE - self.getMatrix()[randSpin[0]][randSpin[1]] * self.getMatrix()[randSpin[0]][randSpin[1]-self.getSideLength()+1]
        randSpinE = randSpinE * (2.0) * self.getCoupling()
            #randSpinE = randSpinE * self.getCoupling()
        return randSpinE
    
    def calcDifferenceInMagnetization(self,randSpin):
        aveSpin = 2*self.getMatrix()[randSpin[0]][randSpin[1]] / math.pow(self.getSideLength(),2)
        return aveSpin

    def accumulateStatistics(self,iStep,energyOfStep,magnetizationOfStep,spinPairCorrOfStep):
        #increment histograms for energy, call sum of spins, etc.
        aveSpinOfStep = [iStep,magnetizationOfStep]
        self.aveSpin.append(aveSpinOfStep)
        latticeEnergyOfStep = [iStep,energyOfStep]
        self.latticeEnergy.append(latticeEnergyOfStep)
        spinPairCorrValue = [iStep,spinPairCorrOfStep]
        self.spinPairCorr.append(spinPairCorrValue)

def MonteCarloSim(nCycles,sideLength,temperature,coupling):
    lattice = Lattice(sideLength)
    lattice.setTemperature(temperature)
    lattice.setCoupling(coupling)
    nSpins = 1 
    acceptedMoves = 0.0
    attemptedMoves = 0.0
    oldEnergy = lattice.calcTotalEnergy()
    oldMagnetization = lattice.calcMagnetization()
    randomSpin = []
    sampleTime = 1000
    distance = 8
    #lattice.printMatrix()
    for i in range(0,int(nCycles)):
        del randomSpin[:]
        deltaE = 0.0
        deltaS = 0.0
        if (i % sampleTime) == 0:
            spinPairCorrOfStep = lattice.spinPair(distance)
            lattice.accumulateStatistics(i,oldEnergy,oldMagnetization,spinPairCorrOfStep)
        for j in range(0,nSpins):
            attemptedMoves = attemptedMoves + 1
            randomSpin.append(lattice.getRandomSpin())
        for j in range(0, len(randomSpin)):   
            lattice.flipSpin(randomSpin[j])
            deltaE = deltaE + lattice.calcDifferenceInEnergy(randomSpin[j])
            deltaS = deltaS + lattice.calcDifferenceInMagnetization(randomSpin[j])
        newEnergy = oldEnergy + deltaE
        newMagnetization = oldMagnetization + deltaS

        if newEnergy > oldEnergy:
            ranNum = random.uniform(0,1)
            if ranNum <= math.exp(-(1.0/lattice.getTemperature())*(newEnergy - oldEnergy)):
                acceptedMoves = acceptedMoves + 1*nSpins
                oldEnergy = newEnergy
                oldMagnetization = newMagnetization
                
            else:
                #newEnergy is to large to be accepted, revert back to oldLattice
                for j in range(0, len(randomSpin)):
                    lattice.flipSpin(randomSpin[j])
                oldEnergy = newEnergy - deltaE
                oldMagnetization = newMagnetization - deltaS
        else:
            #newEnergy is less than oldEnergy therefore keep the new configuration
            acceptedMoves = acceptedMoves + 1*nSpins
            oldEnergy = newEnergy
            oldMagnetization = newMagnetization
            
        if (i % 100000) ==0:
            print "nSpins used = ",nSpins
            print "efficiency = ", acceptedMoves / attemptedMoves
            if (acceptedMoves / attemptedMoves) > 0.50:
                nSpins = nSpins + 1
            else:
                if nSpins > 1:
                    nSpins = nSpins - 1
            print "energy of lattice = ", oldEnergy
            print "magnetization of lattice = ", oldMagnetization
            print "%complete = ", (i/nCycles)*100
            acceptedMoves = 0.0
            attemptedMoves = 0.0
         
        if (i % (nCycles-1)) == 0 and (i != 0):
            #lattice.printMatrix()
            spinPairCorrFile = open('spinPairCorr','w')
            aveSpinFile = open('aveSpin', 'w')
            latticeEnergyFile = open('latticeEnergy', 'w')
            for i in range(0, len(lattice.latticeEnergy)):
                aveSpinFile.write(str(lattice.aveSpin[i][0])+"\t"+str(lattice.aveSpin[i][1]))
                aveSpinFile.write("\n")
                latticeEnergyFile.write(str(lattice.latticeEnergy[i][0])+"\t"+str(lattice.latticeEnergy[i][1]))
                latticeEnergyFile.write("\n")
                spinPairCorrFile.write(str(lattice.spinPairCorr[i][0])+"\t"+str(lattice.spinPairCorr[i][1]))
                spinPairCorrFile.write("\n")
                                       
            #print lattice.aveSpin

       
MonteCarloSim(1E6,20,2.2,1)


def TestRanNumGen(nTrials, nBins):
    #Create and zero histrogram for random numbers
    ranNumHist = []
    for i in range(0, nBins):
        ranNumHist.append(0)
    #Populate histogram by generating random numbers and casting into histogram
    for i in range(0, nTrials):
        ranNumHist[int(random.uniform(0,1)*nBins)]+=1
    #Plot results
    binIndex = []
    binCount = []
    for i in range(0, len(ranNumHist)):
        binIndex.append(i*(1/float(nBins)))
        binCount.append((ranNumHist[i]/float(nTrials)*nBins))
    plt.plot(binIndex,binCount)
    plt.ylabel('bin weighted fractional selection')
    plt.xlabel('x-value')
    plt.show()
    return 


#TestRanNumGen(1000000,100)
