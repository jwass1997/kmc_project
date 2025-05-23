#include <algorithm>
#include <cmath>

#include "Random.h"
#include "State.h"
#include "FEMmethods.h"
#include "Configuration.h"

State::State()
    : nAcceptors(200)
    , nDonors(3)
    , nElectrodes(8)
    , numOfSites(nAcceptors + nElectrodes)
    , radius(150.0)
    , nu0(1.0)
    , a(20.0)
    , T(77.0)
    , energyDisorder(0.05*e / kb*T)
    , minHopDistance(3.0)
    , maxHopDistance(60.0)
    , noDimension(true)
{
    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    occupationOfStates.resize(nAcceptors, 0);

    initRandomState();
}

State::State(Configuration& config) {

    initStateFromConfig(config);
}

void State::initRandomState() {

    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

    if (noDimension) {
        radius = radius / R;
    }

    for (int i = 0; i < nAcceptors; ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        acceptorCoordinates[i*2] = randomR*std::cos(randomPhi);
        acceptorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }

    for (int i = 0; i < nDonors;  ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        donorCoordinates[i*2] = randomR*std::cos(randomPhi);
        donorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }  
}

void State::initStateFromConfig(Configuration& config) {

    nAcceptors = config.nAcceptors;
    nDonors = config.nDonors;
    nElectrodes = config.nElectrodes;
    numOfSites = config.numOfSites;
    radius = config.radius;
    nu0 = config.nu0;
    a = config.a;
    T = config.T;
    energyDisorder = config.energyDisorder;
    R = config.R;
    A0 = config.A0;
    electrodeWidth = config.electrodeWidth;
    minHopDistance = config.minHopDistance;
    maxHopDistance = config.maxHopDistance;

    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    occupationOfStates.resize(nAcceptors, 0);
}

void State::initContainers() {

    std::vector<std::array<double,2>> allCoordinates(numOfSites);

    for (int i = 0; i < nAcceptors; ++i) {
        allCoordinates[i][0] = acceptorCoordinates[i*2];
        allCoordinates[i][1] = acceptorCoordinates[i*2 + 1];
    }

    for (int i = 0; i < nElectrodes; ++i) {
        allCoordinates[i+nAcceptors][0] = electrodeCoordinates[i*2];
        allCoordinates[i+nAcceptors][1] = electrodeCoordinates[i*2 + 1];
    }

    numOfNeighbours.resize(numOfSites);
    int totalNumOfEvents = 0;
    for (int i = 0; i < numOfSites; ++i) {
        distanceMatrix[i*numOfSites + i] = 0.0;
        for (int j = i + 1; j < numOfSites; ++j) {
            double Dx = allCoordinates[i][0] - allCoordinates[j][0];
            double Dy = allCoordinates[i][1] - allCoordinates[j][1];
            double distance = std::sqrt(Dx*Dx + Dy*Dy);
            distanceMatrix[i*numOfSites + j] = distance;
            distanceMatrix[j*numOfSites + i] = distance;
            if ((distance > minHopDistance) && (distance < maxHopDistance)) {
                totalNumOfEvents++;
                numOfNeighbours[i]+=1;
                numOfNeighbours[j]+=1;
            }
        }  
    }

    jaggedArrayLengths.resize(numOfSites+1);
    jaggedArrayLengths[0] = 0;
    for (int i = 0; i < numOfSites; ++i) {
        jaggedArrayLengths[i+1] = jaggedArrayLengths[i] + numOfNeighbours[i];
    }

    std::vector<int> writePtr(numOfSites);
    for (int i = 0; i < numOfSites; ++i) {
        writePtr[i] = jaggedArrayLengths[i];
    }

    constantTransitionRates.resize(2*totalNumOfEvents);
    dynamicalTransitionRates.resize(2*totalNumOfEvents);
    neighbourIndices.resize(2*totalNumOfEvents);

    for (int i = 0; i < numOfSites; ++i) {
        for (int j = i+1; j < numOfSites; ++j) {
            double distance =  distanceMatrix[i*numOfSites + j];
            if (distance > minHopDistance && distance < maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                neighbourIndices[indexIJ] = j;
                neighbourIndices[indexJI] = i;
                constantTransitionRates[indexIJ] = nu0*fastExp(-2.0*distance / a);
                constantTransitionRates[indexJI] = nu0*fastExp(-2.0*distance / a);
            }
        }
    }

    for (int i = 0; i < nAcceptors; ++i) {
        inverseAcceptorDistances[i*nAcceptors + i] = 0.0;
        for (int j = i+1; j < nAcceptors; ++j) {
            double inverseDistance = 1.0 / distanceMatrix[i*numOfSites + j];
            inverseAcceptorDistances[i*nAcceptors + j] = inverseDistance;
            inverseAcceptorDistances[j*nAcceptors + i] = inverseDistance;
        }
    } 
}

void State::initSiteEnergies(FiniteElementeCircle& femSolver) {

    std::vector<double> inverseDistances(nAcceptors, 0.0);
    // Acc-Don interaction + random energy + potential energy (for acceptors only)
    for(int i = 0; i < nAcceptors; ++i) {
        double potentialEnergy = femSolver.getPotential(
            acceptorCoordinates[i*2], 
            acceptorCoordinates[i*2 + 1]
        )*e / kb*T;
        //std::cout << potentialEnergy << "\n";
		double sumOfInverseDistances = 0.0;
		for(int j = 0; j <  nDonors; j++) {
			sumOfInverseDistances += 1.0 / calculateDistance(
                acceptorCoordinates[i*2], 
                donorCoordinates[j*2], 
                acceptorCoordinates[i*2 + 1], 
                donorCoordinates[j*2 + 1]
            );
		}

		inverseDistances[i] = sumOfInverseDistances;
		acceptorDonorInteraction[i] = A0*sumOfInverseDistances;

		if(energyDisorder != 0.0) {
			double randomEnergy = normalDist(0.0, energyDisorder);	
			randomEnergies[i] = randomEnergy;		
		}
        stateEnergies[i] = potentialEnergy + acceptorDonorInteraction[i] + randomEnergies[i];
	}
    // Potential energy (for electrodes only)
	for(int i = 0; i < nElectrodes; ++i) {
		stateEnergies[i] = femSolver.getPotential(electrodeCoordinates[i*2], electrodeCoordinates[i*2 + 1])*e / kb*T;
	}
    // Acc-Acc interaction
    for (int i = 0; i < nAcceptors; ++i) {
        for (int j = 0; j < nAcceptors; ++j) {
            if (i != j) {
                acceptorInteraction[i] += (1 - occupationOfStates[j]) * inverseAcceptorDistances[i*nAcceptors + j];
            }
        }
        stateEnergies[i] += - A0*acceptorInteraction[i];
    }
}

void State::initOccupiedSites() {
    if (nDonors >= nAcceptors) {
        throw std::invalid_argument("initOccupiedStates: Number of acceptors can not be equal or smaller than number of donors!");
    }

    std::vector<int> randomVector(nAcceptors, 0);

    for(int i = 0; i < nAcceptors; ++i) {
        randomVector[i] = i;
    }
    std::shuffle(randomVector.begin(), randomVector.end(), rng);
    for (int i = 0; i < nAcceptors - nDonors; ++i) {
        occupationOfStates[randomVector[i]] = 1;
    }
}

void State::initOccupiedSites() {
    
    int stateToOccupy = 0;
    for (int i = 0; i < nDonors; ++i) {
        stateToOccupy = randomInt(0, nAcceptors);
        occupationOfStates[stateToOccupy] = 1;
    }
}

void State::initOccupiedSitesFromConfig(Configuration& config) {


}