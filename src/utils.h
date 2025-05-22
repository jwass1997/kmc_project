#pragma once

#include <iostream>
#include <random>
#include <cstdint>
#include <cstring>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>

#include "cnpy.h"

class SystemGraph;
class Simulator;

struct Electrode {
    double angularPosition;
    double voltage;
};

extern std::mt19937 rng;

inline constexpr double kb = 1.380649e-23;

inline constexpr double e = 1.602176634e-19;

inline constexpr double PI = 3.1415926535897;

inline constexpr double eps0 = 8.854187817620389e-12;

inline constexpr double epsr = 10.0;

inline constexpr double kbT = 1.0630997e-21;

inline double fastExp(double x);

int argParser(int argc, char* argv[]);

double sampleFromUniformDistribution(double min, double max);

double sampleFromNormalDistribution(double mean, double standardDeviation);

double calculateDistance(double coordinateX1, double coordinateX2, double coordinateY1, double coordinateY2);

void createDirectoryFromStringPath(const std::string& path, const std::string& directoryName);

void recordDevice(const std::string& ID, int equilibriumSteps, int numOfSteps, const std::string& defaultConfigs, const std::string& saveFolderPath);

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

inline std::vector<double> sampleVoltageSetting(
    int numOfElectrodes,
    double minVoltage,
    double maxVoltage) {

        std::vector<double> voltageSetting(numOfElectrodes, 0.0);

        for (auto& voltage : voltageSetting) {
            voltage = sampleFromUniformDistribution(minVoltage, maxVoltage);
        }
        
        return voltageSetting;
};

double IVPoint(
    std::vector<double> voltageSetting,
    int numOfDevices,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfig
);

double currentFromVoltageCombination(
    Simulator& simulator,
    std::vector<double> voltageSetting,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfig
);

void createBatchOfSingleSystem(
    int batchSize, 
    int outputElectrodeIndex,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps, 
    int simulationSteps, 
    int numOfIntervals,
    const std::string& defaultConfigs, 
    const std::string& saveFolder, 
    const std::string& batchID
);