#include <algorithm>
#include <cmath>

#include "State.h"
#include "FEMmethods.h"
#include "Configuration.h"

State::State() {

    numOfSites = nAcceptors + nElectrodes;

    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    occupationOfStates.resize(nAcceptors, 0);
}

State::State(const std::string& configPath) {



}

void State::initRandomState() {

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

void State::initStateFromConfig(Configuration config) {
    
}