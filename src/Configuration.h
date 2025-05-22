#pragma once

#include <filesystem>

class Configuration {
    public:

        Configuration();

        Configuration(const std::string& configPath);

        int nAcceptors;
        int nDonors;
        int nElectrodes;
        
        int numOfSites;

        double radius;

        double nu0;
        double a;
        double T;

        double energyDisorder;
        double R;
        double A0;

        double electrodeWidth;

        double minHopDistance;
        double maxHopDistance;

        std::vector<double> acceptorCoords;

        std::vector<double> donorCoords;

        std::vector<double> electrodeCoords;

        std::vector<int> siteOccupation;
    private:

        std::filesystem::path getConfigFilePath(const std::string& folder, const std::string& file);
};