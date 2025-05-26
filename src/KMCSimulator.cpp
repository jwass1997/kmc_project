#include <vector>
#include <cmath>

#include "KMCParameters.h"
#include "KMCSimulator.h"
#include "Random.h"
#include "State.h"
#include "utils.h"

struct Results {
    
};

KMCSimulator::KMCSimulator()
{
    std::cout << "KMCSimulator(): Empty constructor should not be used" << "\n";
}

KMCSimulator::KMCSimulator(State& state)
{
    numOfNeighbors = state.numOfNeighbours;
    jaggedArrayLengths = state.jaggedArrayLengths;
    neighborIndices = state.neighbourIndices;
    totalNumOfEvents = state.totalNumOfEvents;
    lastHopIndices.resize(2, 0);

    constantTransitionRates.resize(2*totalNumOfEvents, 0.0);
    dynamicalTransitionRates.resize(2*totalNumOfEvents, 0.0);
    aggregatedTransitionRates.resize(2*totalNumOfEvents, 0.0);

    std::vector<int> writePtr(state.numOfSites);
    for (int i = 0; i < state.numOfSites; ++i) {
        writePtr[i] = jaggedArrayLengths[i];
    }

    for (int i = 0; i < state.numOfSites; ++i) {
        for (int j = i+1; j < state.numOfSites; ++j) {
        double distance =  state.distanceMatrix[i*state.numOfSites + j];
            if (distance > state.minHopDistance && distance < state.maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                constantTransitionRates[indexIJ] = state.nu0*fastExp(-2.0*distance / state.a);
                constantTransitionRates[indexJI] = state.nu0*fastExp(-2.0*distance / state.a);
            }
        }
    }
}

void KMCSimulator::initKMCSimulator(State& state) {

    numOfNeighbors = state.numOfNeighbours;
    jaggedArrayLengths = state.jaggedArrayLengths;
    neighborIndices = state.neighbourIndices;
    totalNumOfEvents = state.totalNumOfEvents;
    lastHopIndices.resize(2, 0);

    constantTransitionRates.resize(2*totalNumOfEvents, 0.0);
    dynamicalTransitionRates.resize(2*totalNumOfEvents, 0.0);
    aggregatedTransitionRates.resize(2*totalNumOfEvents, 0.0);

    std::vector<int> writePtr(state.numOfSites);
    for (int i = 0; i < state.numOfSites; ++i) {
        writePtr[i] = jaggedArrayLengths[i];
    }

    for (int i = 0; i < state.numOfSites; ++i) {
        for (int j = i+1; j < state.numOfSites; ++j) {
        double distance =  state.distanceMatrix[i*state.numOfSites + j];
            if (distance > state.minHopDistance && distance < state.maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                constantTransitionRates[indexIJ] = state.nu0*fastExp(-2.0*distance / state.a);
                constantTransitionRates[indexJI] = state.nu0*fastExp(-2.0*distance / state.a);
            }
        }
    }
}

void KMCSimulator::updateTransitionRates(State& state) {

    totalSumOfRates = 0.0;

    for (int i = 0; i < state.numOfSites; ++i) {
		for (int k = jaggedArrayLengths[i]; k < jaggedArrayLengths[i+1]; ++k) {
            int partner = neighborIndices[k];
			// Electrode - Acceptor
			if (i >= state.nAcceptors && partner < state.nAcceptors) {
				if(state.currentOccupation[partner] == 1) {
					dynamicalTransitionRates[k] = 0.0;
				}
				else {
					double deltaE = state.siteEnergies[partner] - state.siteEnergies[i];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = state.nu0;
					} 
					else {
						dynamicalTransitionRates[k] = state.nu0*fastExp(-deltaE);
					} 
				}
			}
			// Acceptor - Electrode
			else if (i < state.nAcceptors && partner >= state.nAcceptors) {
				if (state.currentOccupation[i] == 0) {
					dynamicalTransitionRates[k] = 0.0;
				} 
				else {
					double deltaE = state.siteEnergies[partner] - state.siteEnergies[i];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = state.nu0;
					}
					else {
						dynamicalTransitionRates[k] = state.nu0*fastExp(-deltaE);
					}
				}
			}
			// Acceptor - Acceptor
			else if (i < state.nAcceptors && partner < state.nAcceptors) {
				if ((state.currentOccupation[i] == 1) && (state.currentOccupation[partner] == 0)) {
					double deltaE = state.siteEnergies[partner] - state.siteEnergies[i] - state.A0 / state.distanceMatrix[i*state.numOfSites + partner];
					if (deltaE < 0.0) {
						dynamicalTransitionRates[k] = state.nu0;
					}
					else {
						dynamicalTransitionRates[k] = state.nu0*fastExp(-deltaE);
					} 
				}
				else {
					dynamicalTransitionRates[k] = 0.0;
				}					
			}
            aggregatedTransitionRates[k] = dynamicalTransitionRates[k]*constantTransitionRates[k];
            totalSumOfRates += aggregatedTransitionRates[k];
		}
	}
}

void KMCSimulator::mcStep(State& state, bool writeData) {

    state.updateSiteEnergies(lastHopIndices);
    updateTransitionRates(state);
    sampleEvent(state);
    state.updateSiteOccupation(lastHopIndices);
    state.increaseStateTime(totalSumOfRates);

    if (writeData) {
        state.eventCounter[lastHopIndices[0]*state.numOfSites + lastHopIndices[1]] += 1;
    }
}

void KMCSimulator::simulate(State& state, int steps, bool reset, bool writeData) {

    if (reset) {
        state.resetState();
    }

    for (int i = 0; i < steps; ++i) {
        mcStep(state, writeData);
    }
}

void KMCSimulator::sampleEvent(State& state) {
    
    double _r = totalSumOfRates*randomDouble01();
    cumulativeSumOfRates = 0.0;
    for (int i = 0; i < state.numOfSites; ++i) {
        int L = jaggedArrayLengths[i];
        int R = jaggedArrayLengths[i+1];
        for (int k = L; k < R; ++k) {
            double rate = aggregatedTransitionRates[k];
            cumulativeSumOfRates+=rate;
            if (cumulativeSumOfRates >= _r) {
                int j = neighborIndices[k];
                lastHopIndices[0] = i;
                lastHopIndices[1] = j;
                return;
            }
        }
    }
}