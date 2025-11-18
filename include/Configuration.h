#pragma once

#include <filesystem>
#include <random>
#include <vector>

#include "utils.h"

bool subset_connected(
    std::vector<double> &subset_coords, 
    std::vector<double> &coords, 
    double min_dist, 
    double max_dist
);

std::vector<double> gm_rej_sampling(
    int num_samples, 
    double radius,
    double min_distance,
    double max_distance,
    GaussianMixture2D& gaussian_mixture,
    std::vector<double>& subset_nodes,
    const int max_attempts
);

struct ConfigurationParams {

    int nAcceptors = 200;
    int nDonors = 3;
    int nElectrodes = 8;

    double radius = 150.0;
    double nu0 = 1.0;
    double a = 20.0;
    double T = 77.0;

    double energyDisorder = 0.01;

    double electrodeWidth = 60.0;
    double minHopDistance = 1.5;
    double maxHopDistance = 60.0;

    bool noDimension = true;

    int Nr = 257;
    int Nt = 1440;

    std::vector<Electrode> electrodeData = {
        {0.0, 0.0},
        {45.0, 0.0},
        {90.0, 0.0},
        {135.0, 0.0},
        {180.0, 0.0},
        {225.0, 0.0},
        {270.0, 0.0},
        {315.0, 0.0}
    };

    std::string distType = "uniform";

    int n_comps;

    std::vector<double> var11;
    std::vector<double> var22;
    std::vector<double> var12;
    std::vector<double> m1;
    std::vector<double> m2;
    std::vector<double> normalized_weights;
};

class Configuration {

    public:

        Configuration() = delete;

        Configuration(
            ConfigurationParams& params,
            uint64_t seed      
        );

        Configuration(
            const std::string& cfg, 
            const std::string& acceptorCfg,
            const std::string& donorCfg,
            const std::string& electrodeCfg
        );

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

        bool noDimension;

        int Nr;
        int Nt;

        std::string distType;

        std::vector<int> occupiedSites;

        std::vector<double> acceptorCoords;

        std::vector<double> donorCoords;

        std::vector<double> electrodeCoords;

        std::vector<int> siteOccupation;

        std::vector<Electrode> electrodeData; 

        double uniform01();

        std::mt19937_64 _rng;
};