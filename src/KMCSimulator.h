#pragma once

#include <vector>
#include <random>
#include <iostream>

class State;
class KMCParameters;

class KMCSimulator {

    public:

        KMCSimulator();

        KMCSimulator(State& state, KMCParameters& kmcParams);

        int numOfSteps;

        int equilibriumSteps;

        void initKMCSimulator(State& state, KMCParameters& kmcParams);
        
        void updateTransitionRates(State& state);

        void sampleEvent(State& state);

        void simulate(State& state, bool reset);

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<double> aggregatedTransitionRates;

        int totalSumOfRates;

        int cumulativeSumOfRates;

        std::vector<int> numOfNeighbors;

        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighborIndices;

        int totalNumOfEvents;

        std::vector<int> lastHopIndices;
};