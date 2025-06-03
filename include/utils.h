#pragma once

#include <iostream>
#include <random>
#include <cstdint>
#include <cstring>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>

#include "cnpy.h"

class KMCSimulator;
class State;

struct Electrode {
    double angularPosition;
    double voltage;
};

inline constexpr double kb = 1.380649e-23;

inline constexpr double e = 1.602176634e-19;

inline constexpr double PI = 3.1415926535897;

inline constexpr double eps0 = 8.854187817620389e-12;

inline constexpr double epsr = 10.0;

inline double fastExp(double x);

inline double fastExp(double x) {
    /**
     * Fast exponential from: https://gist.github.com/jrade/293a73f89dfef51da6522428c857802d
     */
    constexpr double a = (1ll << 52) / 0.6931471805599453;
    constexpr double b = (1ll << 52) * (1023 - 0.04367744890362246);
    x = a * x + b;

    // Remove these lines if bounds checking is not needed
    constexpr double c = (1ll << 52);
    constexpr double d = (1ll << 52) * 2047;
    if (x < c || x > d)
        x = (x < c) ? 0.0 : d;

    // With C++20 one can use std::bit_cast instead
    uint64_t n = static_cast<uint64_t>(x);
    memcpy(&x, &n, 8);
    return x;
};

double calculateDistance(double coordinateX1, double coordinateX2, double coordinateY1, double coordinateY2);

void createDirectoryFromStringPath(const std::string& path, const std::string& directoryName);

void singleRun(
    const std::string& ID, 
    int equilibriumSteps, 
    int numOfSteps, 
    const std::string& defaultConfigs, 
    const std::string& saveFolderPath
);

double calculateCurrent(
    State& state,
    KMCSimulator& kmc,
    int electrodeIdx,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals
);

void createDatapoint(
    const std::string& name,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& saveFolderPath
);

void singleStateBatch(
    int batchSize,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& configs,
    const std::string& save,
    const std::string& batchName
);

int argParser(int argc, char* argv[]);