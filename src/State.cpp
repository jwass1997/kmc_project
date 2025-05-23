#include <algorithm>
#include <cmath>

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