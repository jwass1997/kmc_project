#pragma once

#include <filesystem>

#include "utils.h"

class Configuration {

    public:

        Configuration();

        Configuration(const std::string& configPath, bool randomCoordinates);

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

        bool noDimension = true;

        std::vector<int> occupiedSites;

        std::vector<double> acceptorCoords;

        std::vector<double> donorCoords;

        std::vector<double> electrodeCoords;

        std::vector<int> siteOccupation;

        std::vector<Electrode> electrodeData; 

    private:
    
        std::filesystem::path getConfigFilePath(const std::string& folder, const std::string& file);
};