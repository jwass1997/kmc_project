#pragma once

#include <vector>
#include <string>
#include <cmath>

#include "utils.h"

class Configuration;
class FiniteElementCircle;

class State {

    public:

        State();

        State(Configuration& config);

        void initRandomState();
        
        void initStateFromConfig(Configuration& config);

        void initContainers();

        void initSiteEnergies(FiniteElementeCircle& femSolver);

        void initOccupiedSites();

        void initOccupiedSitesFromConfig(Configuration& config);

        int nAcceptors;

        int nDonors;

        int nElectrodes;

        int numOfSites;

        double radius;

        double nu0;
        double a;
        double T;

        double energyDisorder;
        double R;
        double A0;

        double electrodeWidth;

        double minHopDistance;
        double maxHopDistance;

        bool noDimension = true;

        std::vector<double> acceptorCoordinates;

        std::vector<double> donorCoordinates;

        std::vector<double> electrodeCoordinates;

        std::vector<double> distanceMatrix;

        std::vector<double> inverseAcceptorDistances;

        std::vector<Electrode* > electrodeData;

        std::vector<int> occupationOfStates;

        std::vector<double> randomEnergies;

        std::vector<double> acceptorDonorInteraction;

        std::vector<double> acceptorInteraction;

        std::vector<double> stateEnergies;        

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<int> numOfNeighbours;
        
        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighbourIndices;
};