#pragma once

#include <vector>
#include <string>
#include <cmath>

#include "utils.h"

class Configuration;
class FiniteElementeCircle;

class State {

    public:

        State();

        State(Configuration& config, FiniteElementeCircle& fem);

        void initRandomState();
        
        void initStateFromConfig(Configuration& config);

        void initContainers();

        void initSiteEnergies(FiniteElementeCircle& femSolver);

        void initOccupiedSites();

        void initOccupiedSitesFromConfig(Configuration& config);

        void updateSiteEnergies(std::vector<int> lastHopIndices);

        void updateSiteOccupation(std::vector<int> lastHopIndices);

        void updateBoundaries(std::vector<double> boundaryValues, FiniteElementeCircle& fem);

        void increaseStateTime(double rate);

        void resetState();

        void updatePotential(FiniteElementeCircle& fem);

        void resetPotential(FiniteElementeCircle& fem);

        int nAcceptors;

        int nDonors;

        int nElectrodes;

        int numOfSites;

        double radius;

        double nu0;
        double a;
        double T;

        double kbT;

        double energyDisorder;
        double R;
        double A0;

        double electrodeWidth;

        double minHopDistance;
        double maxHopDistance;

        bool noDimension = true;

        double stateTime;

        std::vector<double> acceptorCoordinates;

        std::vector<double> donorCoordinates;

        std::vector<double> electrodeCoordinates;

        std::vector<double> distanceMatrix;

        std::vector<double> inverseAcceptorDistances;

        std::vector<Electrode> electrodeData;

        std::vector<int> initialOccupation;

        std::vector<int> currentOccupation;

        std::vector<double> randomEnergies;

        std::vector<double> acceptorDonorInteraction;

        std::vector<double> acceptorInteraction;

        std::vector<double> initialSiteEnergies;  
        
        std::vector<double> initialPotential;
        
        std::vector<double> currentPotential;

        std::vector<double> siteEnergies;

        std::vector<int> numOfNeighbours;
        
        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighbourIndices;

        int totalNumOfEvents;

        std::vector<int> eventCounter;
};