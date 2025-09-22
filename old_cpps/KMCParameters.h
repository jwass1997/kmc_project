#pragma once

#include <string>
#include <filesystem>

class KMCParameters {
    
    public:

        KMCParameters();
        KMCParameters(const std::string& paramsPath);

        int numOfSteps;

        int equilibriumSteps;

    private:

        std::filesystem::path getConfigFilePath(const std::string& folder, const std::string& file);
};