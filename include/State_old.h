#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <random>

#include "utils.h"
#include "mfem.hpp"

class Configuration;
class FiniteElementeCircle;

class State {

    public:

        State();

        State(State const& other);

        State(Configuration& config, uint64_t seed);

        void initContainers();

        void initPotential();

        void initSiteEnergies();

        void initOccupiedSites();

        void updateSiteEnergies(std::vector<int> lastHopIndices);

        void updateSiteOccupation(std::vector<int> lastHopIndices);

        void updateBoundaries(std::vector<double> boundaryValues);

        void updateBoundaries(mfem::GridFunction const& solutionVector);

        mfem::GridFunction getSolutionVector() const;

        void increaseStateTime(double rate);

        void resetEventCounter();

        void resetState();

        double uniform01();

        double normal(double mu, double sigma);

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

        int femRes;

        std::unique_ptr<FiniteElementeCircle> femSolver;

        std::mt19937_64 _rng;
};