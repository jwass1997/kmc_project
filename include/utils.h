#pragma once

#include <array>
#include <cwchar>
#include <linux/limits.h>
#include <random>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <cmath>

#include "cnpy.h"

class KMCSimulator;
class State;

struct Electrode {
    double angularPosition;
    double voltage;
};

struct Gaussian2D {
    double mx, my; 
    double sxy, sxx, syy;

    double det;
    double inv00, inv01, inv10, inv11;
    double norm_const;

    Gaussian2D(double mx, double my, double sxy, double sxx, double syy) 
    : mx(mx)
    , my(my)
    , sxy(sxy)
    , sxx(sxx)
    , syy(syy)
    {
        det = sxx * syy - sxy * sxy;

        if (det <= 0.0) 
        {
            throw std::runtime_error("Covariance is not not positive definite.");
        }
            
        inv00 = syy / det;
        inv01 = -sxy  / det;
        inv10 = -sxy / det;
        inv11 = sxx / det;

        norm_const = 1.0 / (2.0*M_PI * std::sqrt(det));
    }

    std::array<double, 2> sample (std::mt19937 &rng) const{
        thread_local std::normal_distribution<double> normal_dist_1D(0.0, 1.0);

        double a = sxx;
        double b = sxy;
        double c = syy;

        double l11 = std::sqrt(a);

        if (l11 <= 0.0)
        {
            throw std::runtime_error("Cholesky failed: non-SPD covariance");
        }
        
        double l21 = b / l11;
        double t = c - l21 * l21;

        if (t <= 0.0)
        {
            throw std::runtime_error("Cholesky failed: non-SPD covariance.");
        }
        
        double l22 = std::sqrt(t);

        double z1 = normal_dist_1D(rng);
        double z2 = normal_dist_1D(rng);

        double x = mx + l11*z1;
        double y = my + l21*z1 + l22*z2;

        return {x, y};
    }
};

struct GaussianMixture2D {

    std::vector<Gaussian2D> components;
    std::vector<double> cum_weights;
    std::vector<double> normalized_weights;

    GaussianMixture2D (const std::vector<Gaussian2D> &components, const std::vector<double> &weights)
    :
    components(components)
    {
        if (components.empty() || weights.size() != components.size())
        {
            throw std::runtime_error("Mismatch between components and weight sizes or empty components.");
        }

        double sum = 0.0;
        for (auto w : weights)
        {
            sum += w;
        }

        if (sum <= 0)
        {
            throw std::runtime_error("Sum of weights must be positive.");
        }

        cum_weights.reserve(weights.size());
        double w_sum = 0.0;
        for (auto w : weights) 
        {
            w_sum += w / sum;
            cum_weights.push_back(w_sum);
            normalized_weights.push_back(w / sum);
        }

        cum_weights.back() = 1.0;
    }

    static std::mt19937 &thread_rng () {
        thread_local std::mt19937 rng{std::random_device{}()};
        return rng;
    }

    std::array<double, 2> sample () const {
        auto &rng = thread_rng();

        std::uniform_real_distribution<double> uni_01(0.0, 1.0);
        double u = uni_01(rng);

        int k = 0;
        while (k + 1 < cum_weights.size() && u > cum_weights[k])
        {
            ++k;
        }

        return components[k].sample(rng);
    }
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

uint64_t mix(uint64_t x);

double calculateDistance(
    double coordinateX1, 
    double coordinateX2, 
    double coordinateY1,
    double coordinateY2
);

double calculateCurrent(
    State& state,
    KMCSimulator& kmc,
    int electrodeIdx,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals
);

void createDirectoryFromStringPath(
    const std::string& path, 
    const std::string& directoryName
);

std::vector<std::vector<double>> scaledLHC(
    std::size_t N,
    std::size_t D,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    unsigned int seed = std::random_device{}()
);

void singleRun(
    const std::string& ID, 
    int eqSteps, 
    int simSteps, 
    int num_intervals,
    std::vector<double> voltages,
    int outputIdx,
    const std::string& cfg, 
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    int seed,
    const std::string& saveFolderPath
);

double singleIVPoint(
    State& initState,
    int outputIdx,
    int numOfTasks,
    int simSteps,
    std::vector<double> voltages
);

void singleIVCurve(
    int numOfPoints,
    int inputIdx,
    int outputIdx,
    double minVoltage,
    double maxVoltage,
    int eqSteps,
    int simSteps,
    int numIntervals,
    int seed,
    std::vector<double> controlVoltages,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolder,
    const std::string& fileName
);

void batchFromSingleState(
    int batchSize,
    double minVoltage,
    double maxVoltage,
    int outputIdx,
    int eqSteps,
    int simSteps,
    int numOfTasks,
    int LHCSeed,
    int threadBaseSeed,
    const std::string& cfg,
    const std::string& accCfg,
    const std::string& donCfg,
    const std::string& eleCfg,
    const std::string& saveFolder,
    const std::string& fileName
);

void batchOfMultipleStates(
    int batchSize,
    double minVoltage, double maxVoltage,
    int nAcceptors, int nElectrodes, int nDonors,
    double radius, double nu0, double a, double T,
    double energyDisorder,
    double electrodeWidth,
    double minHopDistance, double maxHopDistance,
    int Nr, int Nt,
    std::string distType, int n_comps,
    int inputIdx, int outputIdx,
    int eqSteps, int simSteps, int numOfTasks,
    int LHCSeed, int threadBaseSeed,
    const std::string& saveFolder, const std::string& fileName
);

void batchOfSamplesWTemp(
    int batchSize,
    double minVoltage, double maxVoltage,
    double Tmin, double Tmax,
    int nAcceptors, int nElectrodes, int nDonors,
    double radius, double nu0, double a,
    double energyDisorder,
    double electrodeWidth,
    double minHopDistance, double maxHopDistance,
    int Nr, int Nt,
    std::string distType, int n_comps,
    int inputIdx, int outputIdx,
    int eqSteps, int simSteps, int numOfTasks,
    int LHCSeed, int threadBaseSeed,
    const std::string& saveFolder, const std::string& fileName
);

int argParser(int argc, char* argv[]);