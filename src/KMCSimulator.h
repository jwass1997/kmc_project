#pragma once

#include <vector>
#include <random>
#include <iostream>

class State;
class KMCParameters;

class KMCSimulator {

    public:

        KMCSimulator();

        KMCSimulator(State& state);

        int numOfSteps;

        int equilibriumSteps;

        void initKMCSimulator(State& state);
        
        void updateTransitionRates(State& state);

        void sampleEvent(State& state);

        void mcStep(State& state, bool writeData);

        void simulate(State& state, int steps, bool reset, bool writeData);

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<double> aggregatedTransitionRates;

        double totalSumOfRates;

        double cumulativeSumOfRates;

        std::vector<int> numOfNeighbors;

        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighborIndices;

        int totalNumOfEvents;

        std::vector<int> lastHopIndices;
};