#pragma once

#include <vector>
#include <random>
#include <iostream>

class State;
class KMCParameters;

class KMCSimulator {

    public:

        KMCSimulator();

        KMCSimulator(State& state, uint64_t seed);

        void initKMCSimulator(State& state);
        
        void updateTransitionRates(State& state);

        void sampleEvent(State& state);

        void mcStep(State& state, bool writeData);

        void simulate(State& state, int steps, bool reset, bool writeData);

        void resetSimulator();

        double uniform01();

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        std::vector<double> aggregatedTransitionRates;

        double totalSumOfRates;

        double cumulativeSumOfRates;

        std::vector<int> numOfNeighbours;

        std::vector<int> jaggedArrayLengths;

        std::vector<int> neighbourIndices;

        int totalNumOfEvents;

        std::vector<int> lastHopIndices;

        std::mt19937_64 _rng;
};