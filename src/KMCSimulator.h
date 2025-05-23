#pragma once

#include <vector>
#include <random>
#include <iostream>

class State;

class KMCSimulator {

    public:

        KMCSimulator(State& state);

    private:

        std::vector<int> lastHopIndices;

        std::vector<double> constantTransitionRates;

        std::vector<double> dynamicalTransitionRates;

        
};