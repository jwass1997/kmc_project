#include <algorithm>
#include <cmath>
#include <random>

#include "State.h"
#include "FEMmethods.h"
#include "Configuration.h"

State::State() {};

State::State(State const& other)
    : acceptorCoordinates(other.acceptorCoordinates)
    , donorCoordinates(other.donorCoordinates)
    , electrodeCoordinates(other.electrodeCoordinates)
    , nAcceptors(other.nAcceptors)
    , nDonors(other.nDonors)
    , nElectrodes(other.nElectrodes)
    , numOfSites(other.numOfSites)
    , radius(other.radius)
    , nu0(other.nu0)
    , a(other.a)
    , T(other.T)
    , kbT(other.kbT)
    , energyDisorder(other.energyDisorder)
    , R(other.R)
    , A0(other.A0)
    , electrodeWidth(other.electrodeWidth)
    , minHopDistance(other.minHopDistance)
    , maxHopDistance(other.maxHopDistance)
    , electrodeData(other.electrodeData)
    , distanceMatrix(other.distanceMatrix)
    , inverseAcceptorDistances(other.inverseAcceptorDistances)
    , currentOccupation(other.currentOccupation)
    , initialOccupation(other.initialOccupation)
    , randomEnergies(other.randomEnergies)
    , acceptorDonorInteraction(other.acceptorDonorInteraction)
    , acceptorInteraction(other.acceptorInteraction)
    , initialSiteEnergies(other.initialSiteEnergies)
    , initialPotential(other.initialPotential)
    , currentPotential(other.currentPotential)
    , siteEnergies(other.siteEnergies)
    , eventCounter(other.eventCounter)
    , jaggedArrayLengths(other.jaggedArrayLengths)
    , neighbourIndices(other.neighbourIndices)
    , numOfNeighbours(other.numOfNeighbours)
    , totalNumOfEvents(other.totalNumOfEvents)
    , stateTime(other.stateTime)
    , femRes(other.femRes)
    , _rng(other._rng)
{   
    auto femCircle = std::make_unique<FiniteElementeCircle>(
        other.radius,
        other.femRes
    );

    for (int i = 0; i < electrodeData.size(); ++i) {
        femCircle->setElectrode(
            0.0,
            electrodeData[i].angularPosition/360.0 * 2.0*M_PI - 0.5*electrodeWidth / radius,
            electrodeData[i].angularPosition/360.0 * 2.0*M_PI + 0.5*electrodeWidth / radius
        );
    }
    this->femSolver = std::move(femCircle);
    femSolver->initRun();

    *(femSolver->solutionVector) = other.getSolutionVector();
}

State::State(Configuration& config, uint64_t seed)
    : acceptorCoordinates(config.acceptorCoords)
    , donorCoordinates(config.donorCoords)
    , electrodeCoordinates(config.electrodeCoords)
    , nAcceptors(config.nAcceptors)
    , nDonors(config.nDonors)
    , nElectrodes(config.nElectrodes)
    , numOfSites(config.numOfSites)
    , radius(config.radius)
    , nu0(config.nu0)
    , a(config.a)
    , T(config.T)
    , kbT(kb*T)
    , energyDisorder(config.energyDisorder)
    , R(config.R)
    , A0(config.A0)
    , electrodeWidth(config.electrodeWidth)
    , minHopDistance(config.minHopDistance)
    , maxHopDistance(config.maxHopDistance)
    , electrodeData(config.electrodeData)
    , femRes(config.femRes)
    , _rng(seed)
{
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    currentOccupation.resize(nAcceptors, 0);
    initialOccupation.resize(nAcceptors, 0);
    randomEnergies.resize(nAcceptors, 0.0);
    acceptorDonorInteraction.resize(nAcceptors, 0.0);
    acceptorInteraction.resize(nAcceptors*nAcceptors, 0.0);
    initialSiteEnergies.resize(nAcceptors+nElectrodes, 0.0);
    initialPotential.resize(nAcceptors+nElectrodes, 0.0);
    currentPotential.resize(nAcceptors+nElectrodes, 0.0);
    siteEnergies.resize(nAcceptors+nElectrodes, 0.0);   
    eventCounter.resize(numOfSites*numOfSites, 0); 

    initContainers();
    initOccupiedSites();
    initPotential();
    initSiteEnergies();

    stateTime = 0.0;
}

void State::initContainers() {

    std::vector<std::array<double, 2>> allCoordinates(numOfSites);

    for (int i = 0; i < nAcceptors; ++i) {
        allCoordinates[i][0] = acceptorCoordinates[i*2];
        allCoordinates[i][1] = acceptorCoordinates[i*2 + 1];
    }

    for (int i = 0; i < nElectrodes; ++i) {
        allCoordinates[i+nAcceptors][0] = electrodeCoordinates[i*2];
        allCoordinates[i+nAcceptors][1] = electrodeCoordinates[i*2 + 1];
    }

    numOfNeighbours.assign(numOfSites, 0);
    totalNumOfEvents = 0;
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

    neighbourIndices.resize(2*totalNumOfEvents);

    for (int i = 0; i < numOfSites; ++i) {
        for (int j = i+1; j < numOfSites; ++j) {
            double distance =  distanceMatrix[i*numOfSites + j];
            if (distance > minHopDistance && distance < maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                neighbourIndices[indexIJ] = j;
                neighbourIndices[indexJI] = i;
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

void State::initPotential() {

    auto femCircle = std::make_unique<FiniteElementeCircle>(
        this->radius, this->femRes
    );

    for (int i = 0; i < electrodeData.size(); ++i) {
        femCircle->setElectrode(
            0.0,
            electrodeData[i].angularPosition/360.0 * 2.0*M_PI - 0.5*electrodeWidth / radius,
            electrodeData[i].angularPosition/360.0 * 2.0*M_PI + 0.5*electrodeWidth / radius
        );
    }

    this->femSolver = std::move(femCircle);
    femSolver->initRun();

    /* for (int i = 0; i < electrodeData.size(); ++i) {
        femSolver->updateElectrodeVoltage(i, electrodeData[i].voltage);
    }

    femSolver->run(); */
}

void State::initSiteEnergies() {
    
    //std::cout << "Entering InitSiteEnergies" << "\n";
    std::vector<double> inverseDistances(nAcceptors, 0.0);
    // Acc-Don interaction + random energy + potential energy (for acceptors only)
    for(int i = 0; i < nAcceptors; ++i) {
        double potentialEnergy = femSolver->getPotential(
            acceptorCoordinates[i*2], 
            acceptorCoordinates[i*2 + 1]
        )*e / kbT;
        //std::cout << "Potential energy: " << potentialEnergy << "for site " << i << "\n";
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
			double randomEnergy = normal(0.0, energyDisorder);		
            //std::cout << "Random contribution:" << randomEnergy << "\n";
			randomEnergies[i] = randomEnergy;	
            //std::cout << "Random energy: " << randomEnergy << "for site " << i << "\n";
		}
        //std::cout << "Donor interaction " << acceptorDonorInteraction[i] << " for site " << i << "\n";
        initialSiteEnergies[i] = potentialEnergy + acceptorDonorInteraction[i] + randomEnergies[i];
        currentPotential[i] += potentialEnergy;
        initialPotential[i] += potentialEnergy;
        //std::cout << "Current Potential: " << currentPotential[i] << " for site: " << i << "\n";
        //std::cout << "Initial Potential: " << currentPotential[i] << " for site: " << i << "\n";  
	}
    // Potential energy (for electrodes only)
	for(int i = 0; i < nElectrodes; ++i) {
        double potentialEnergy = femSolver->getPotential(electrodeCoordinates[i*2], electrodeCoordinates[i*2 + 1])*e / kbT;
		initialSiteEnergies[i + nAcceptors] += potentialEnergy;
        currentPotential[i + nAcceptors] += potentialEnergy;
        initialPotential[i + nAcceptors] += potentialEnergy;
	}
    // Acc-Acc interaction
    for (int i = 0; i < nAcceptors; ++i) {
        for (int j = 0; j < nAcceptors; ++j) {
            if (i != j) {
                acceptorInteraction[i] += (1 - currentOccupation[j]) * inverseAcceptorDistances[i*nAcceptors + j];
            }
        }
        initialSiteEnergies[i] += - A0*acceptorInteraction[i];
        //std::cout << "Acceptor interaction " << acceptorInteraction[i] << " for site " << i << "\n";
    }
    //std::cout << "Leaving InitSiteEnergies" << "\n";
}

void State::initOccupiedSites() {

    if (nDonors >= nAcceptors) {
        throw std::invalid_argument("[initOccupiedStates]: Number of acceptors can not be equal or smaller than number of donors!");
    }

    std::vector<int> randomVector(nAcceptors, 0);

    for(int i = 0; i < nAcceptors; ++i) {
        randomVector[i] = i;
    }
    std::shuffle(randomVector.begin(), randomVector.end(), _rng);
    std::fill(currentOccupation.begin(), currentOccupation.end(), 0);
    for (int i = 0; i < nAcceptors - nDonors; ++i) {
        currentOccupation[randomVector[i]] = 1;
    }

    initialOccupation = currentOccupation;
}

void State::updateSiteEnergies(std::vector<int> lastHopIndices) {
    
    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] < nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[1]) {
                acceptorInteraction[i] -= 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[1]];
            }
            if (i != lastHopIndices[0]) {
                acceptorInteraction[i] += 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[0]];
            }
        }
    }
    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] >= nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[0]) {
                acceptorInteraction[i] += 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[0]];
            }
        }
    }
    if (lastHopIndices[0] >= nAcceptors && lastHopIndices[1] < nAcceptors) {
        for (int i = 0; i < nAcceptors; ++i) {
            if (i != lastHopIndices[1]) {
                acceptorInteraction[i] -= 1.0*inverseAcceptorDistances[i*nAcceptors + lastHopIndices[1]];
            }
        }
    }

    for (int i = 0; i < nAcceptors; ++i) {
        siteEnergies[i] = currentPotential[i] + acceptorDonorInteraction[i] + randomEnergies[i] - A0*acceptorInteraction[i];
    }
}

void State::updateSiteOccupation(std::vector<int> lastHopIndices) {

    if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] < nAcceptors) {
		currentOccupation[lastHopIndices[0]] = 0;
		currentOccupation[lastHopIndices[1]] = 1;
	}
	if (lastHopIndices[0] < nAcceptors && lastHopIndices[1] >= nAcceptors) {
		currentOccupation[lastHopIndices[0]] = 0;
	}
	if (lastHopIndices[0] >= nAcceptors) {
		if(lastHopIndices[1] < nAcceptors) {
			currentOccupation[lastHopIndices[1]] = 1;
		}
	}
}

void State::updateBoundaries(std::vector<double> boundaryValues) {

    //std::cout << "Entering updateBoundaries" << "\n";
    if (boundaryValues.size() > nElectrodes) {
        throw std::invalid_argument("[updateBoundaries]: Too many boundary values");
    }

    for (int i = 0; i < numOfSites; ++i) {
        siteEnergies[i] -= currentPotential[i];
        //std::cout << "OLD POTENTIAL: " << currentPotential[i] << "\n";
    }

    for (int bdrVal = 0; bdrVal < boundaryValues.size(); ++bdrVal) {
        femSolver->updateElectrodeVoltage(bdrVal, boundaryValues[bdrVal]);
    }

    femSolver->run();
    
    for (int i = 0; i < nAcceptors; ++i) {
        currentPotential[i] = femSolver->getPotential(acceptorCoordinates[i*2], acceptorCoordinates[i*2 + 1])*e / kbT;
        //std::cout << "NEW ACCEPTOR POTENTIAL: " << currentPotential[i] << "\n";
        siteEnergies[i] += currentPotential[i];
    }
    for (int i = 0; i < nElectrodes; ++i) {
        int idx = nAcceptors + i;
        currentPotential[idx] = femSolver->getPotential(electrodeCoordinates[i*2], electrodeCoordinates[i*2 + 1])*e / kbT;
        //std::cout << "NEW ELECTRODE POTENTIAL: " << currentPotential[idx] << "\n";
        siteEnergies[idx] += currentPotential[idx]; 
    }
    //std::cout << "Leaving updateBoundaries" << "\n";
}

void State::updateBoundaries(mfem::GridFunction const& solutionVector) {

}

mfem::GridFunction State::getSolutionVector() const {

    return (*femSolver->solutionVector);
}

void State::increaseStateTime(double rate) {

    double _r = uniform01() + 1e-12;
    double _dt = -(std::log(_r) / rate);

    stateTime += _dt;
}

void State::resetEventCounter() {

    std::fill(eventCounter.begin(), eventCounter.end(), 0);
}

void State::resetState() {

    stateTime = 0.0;

    siteEnergies = initialSiteEnergies;
    currentPotential = initialPotential;
    currentOccupation = initialOccupation;
    std::fill(eventCounter.begin(), eventCounter.end(), 0);
}

double State::uniform01() {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    return U(_rng);
}

double State::normal(double mu, double sigma) {
    std::normal_distribution<double> N(mu, sigma);
    return N(_rng);
}