#pragma once

#include <vector>
#include <random>
#include <iostream>

class State;

class KMCSimulator {

    public:

        KMCSimulator(State& state);

        void initKMCSimulator(State& state);

        void updateTransitionRates(State& state);

        void updateSiteOccupation(State& state);

        void sampleEvent(State& state);

        void simulate(State& state);

        std::vector<double> currentSiteEnergies;     

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