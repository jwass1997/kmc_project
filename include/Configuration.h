#pragma once

#include <filesystem>

#include "utils.h"

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

    int noDimension = 1;

    int femRes = 1e4;

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
    /* Only used if mixed */
    double epsilon = 0.5;
};

class Configuration {

    public:

        Configuration();

        Configuration(
            const ConfigurationParams& params            
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

        int femRes;

        std::string distType;

        double epsilon;

        std::vector<int> occupiedSites;

        std::vector<double> acceptorCoords;

        std::vector<double> donorCoords;

        std::vector<double> electrodeCoords;

        std::vector<int> siteOccupation;

        std::vector<Electrode> electrodeData; 
};