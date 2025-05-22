#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "Configuration.h"

Configuration::Configuration() {

}

Configuration::Configuration(const std::string& configPath) {

    auto config = getConfigFilePath(configPath, "config.txt");
    auto acceptorConfig = getConfigFilePath(configPath, "acceptorConfig.txt");
    auto donorConfig = getConfigFilePath(configPath, "donorConfig.txt");
    auto electrodeConfig = getConfigFilePath(configPath, "electrodeConfig.txt");

    std::ifstream configFile(config);
    std::ifstream acceptorFile(acceptorConfig);
    std::ifstream donorFile(donorConfig);
    std::ifstream electrodeFile(electrodeConfig);

    if (!acceptorFile.is_open()) {
        std::cerr << "No such file: " << acceptorConfig << "\n";
    }
}

std::filesystem::path Configuration::getConfigFilePath(const std::string& folder, const std::string& file) {

    return std::filesystem::path(folder) / file;

}