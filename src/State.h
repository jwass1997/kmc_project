#pragma once

#include <vector>
#include <string>
#include <cmath>

#include "utils.h"

class FiniteElementCircle;

class State {

    public:
        State();

        State(const std::string& configPath);

        std::filesystem::path getConfigFilePath(const std::string& folderPath, const std::string& fileName);

        void initRandomState();
        
        void initStateFromConfig(Configuration config);

        void initContainers();

        void initSiteEnergies();

        void initOccupiedSites();

        int nAcceptors = 200;

        int nDonors = 3;

        int nElectrodes = 8;

        int numOfSites;

        double radius = 150.0;

        double nu0 = 1.0;

        double a = 20.0;

        double T = 77.0;

        double energyDisorder = 0.05*e / kbT;

        double R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

        double A0 = (e*e) / (4.0*kbT*PI*eps0*epsr*1e-9*R);

        double electrodeWidth = 60.0;

        double minHopDistance = 3.0;

        double maxHopDistance = 60.0;

        std::vector<double> acceptorCoordinates;

        std::vector<double> donorCoordinates;

        std::vector<double> electrodeCoordinates;

        std::vector<double> distanceMatrix;

        std::vector<double> inverseAcceptorDistances;

        std::vector<Electrode* > electrodeData;

        std::vector<int> occupationOfStates;

        /* std::vector<double> randomEnergies;

        std::vector<double> acceptorDonorInteraction;

        std::vector<double> acceptorInteraction;

        std::vector<double> stateEnergies;        

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<double> aggregatedTransitionRates; */

        std::vector<int> numOfNeighbours;
        
        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighbourIndices;

        std::vector<int> lastHopIndices;     
};