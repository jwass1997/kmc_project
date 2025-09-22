#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>

#include "KMCParameters.h"

KMCParameters::KMCParameters()
{
    std::cout << "KMCParameters(): Empty constructor should not be used" << "\n";
}

KMCParameters::KMCParameters(const std::string& paramsPath)
{
    auto kmcConfig = getConfigFilePath(paramsPath, "kmc_config.txt");

    std::ifstream kmcFile(kmcConfig);

    if (!kmcFile.is_open()) {
        std::cerr << "KMCParameters: No such file: " << kmcConfig << "\n";
    }

    else {
        std::string line;
        while (getline(kmcFile, line)) {
            if (line.empty() || line[0]=='#') {
                continue;
            }
            else {
                std::stringstream ss(line);
                std::string key, value;

                if (key == "numOfSteps") {
                    numOfSteps = std::stoi(value);
                }
                else if (key == "equilibriumSteps") {
                    equilibriumSteps = std::stoi(value);
                }
            }
        }
        kmcFile.close();
    }
}

std::filesystem::path KMCParameters::getConfigFilePath(const std::string& folder, const std::string& file) {

    return std::filesystem::path(folder) / file;

}