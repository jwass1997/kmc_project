#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <chrono>

#include "Configuration.h"
#include "Random.h"

Configuration::Configuration() 
{
    std::cout << "Configuration(): Empty constructor should not be called!" << "\n";
}

Configuration::Configuration(const std::string& configPath, bool randomCoordinates) {

    auto config = getConfigFilePath(configPath, "config.txt");
    auto acceptorConfig = getConfigFilePath(configPath, "acceptors.txt");
    auto donorConfig = getConfigFilePath(configPath, "donors.txt");
    auto electrodeConfig = getConfigFilePath(configPath, "electrodes.txt");
    /* std::cout
    << "  CWD:    " << std::filesystem::current_path() << "\n"
    << "  Target: " << std::filesystem::absolute(config) << "\n"
    << "  Exists: " << std::boolalpha << std::filesystem::exists(config) << "\n"; */
    std::ifstream configFile(config);
    std::ifstream acceptorFile(acceptorConfig);
    std::ifstream donorFile(donorConfig);
    std::ifstream electrodeFile(electrodeConfig);

    if (!configFile.is_open()) {
        std::cerr << "No such file: " << config << "\n";
    }
    else {
        std::string line;
        while (getline(configFile, line)) {
            if (line.empty() || line[0]=='#') {
                continue;
            }
            else {
                std::stringstream ss(line);
                std::string key, value;

                ss >> key >> value;

                if (key == "nAcceptors") {
                    nAcceptors = std::stoi(value);
                }
                else if (key == "nDonors") {
                    nDonors = std::stoi(value);                   
                }
                else if (key == "nElectrodes") {
                    nElectrodes = std::stoi(value);
                }
                else if (key == "radius") {
                    radius = std::stod(value);
                }
                else if (key == "nu0") {
                    nu0 = std::stod(value);
                }
                else if (key == "a") {
                    a = std::stod(value);
                }
                else if (key == "T") {
                    T = std::stod(value);
                }
                else if (key == "energyDisorder") {
                    energyDisorder = std::stod(value)*e / (kb*T);
                }
                else if (key == "electrodeWidth") {
                    electrodeWidth = std::stod(value);
                }
                else if (key == "minHopDistance") {
                    minHopDistance = std::stod(value);
                }
                else if (key == "maxHopDistance") {
                    maxHopDistance = std::stod(value);
                }
                else if (key == "noDimension") {
                    noDimension = std::stoi(value);
                }
            }
        }
        configFile.close();
    }

    numOfSites = nAcceptors + nElectrodes;
    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));
    kbT = kb*T;
    A0 = (e*e) / (4.0*kb*T*PI*eps0*epsr*1e-9);

    if (noDimension) {
        radius = radius / R;
        a = a / R;
        A0 = A0 / R;
        electrodeWidth = electrodeWidth / R;
        minHopDistance = minHopDistance / R;
        maxHopDistance = maxHopDistance / R;
    }

    // Electrodes
    if (!electrodeFile.is_open()) {
            std::cerr << "No such file: " << electrodeConfig<< "\n";
    }
    
    else {
        std::string line;
        
        while(getline(electrodeFile, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            else {
                std::stringstream ss(line);
                std::string angleStr;
                std::string vStr;

                double angularPosition, voltage;
                ss >> angleStr >> vStr;

                Electrode newElectrode{};

                newElectrode.angularPosition = std::stod(angleStr);
                newElectrode.voltage = std::stod(vStr);

                double phi = (2.0*M_PI*newElectrode.angularPosition) / 360.0;
                double x = radius*std::cos(phi);
                double y = radius*std::sin(phi);
                electrodeCoords.push_back(x);
                electrodeCoords.push_back(y);      
                
                electrodeData.push_back(newElectrode);
            }
        }
        electrodeFile.close();
    }
    // Coordinates
    if (!randomCoordinates) {
        if (!acceptorFile.is_open()) {
            std::cerr << "No such file: " << acceptorConfig << "\n";
        }
        else {
            std::string line;

            while (getline(acceptorFile, line)) {
                if (line.empty() || line[0] == '#') {
                    continue;
                }
                else {
                    std::stringstream ss(line);

                    double coordX, coordY;
                    ss >> coordX >> coordY;
                    
                    if (noDimension) {
                        coordX = coordX / R;
                        coordY = coordY / R;
                    }

                    acceptorCoords.push_back(coordX);
                    acceptorCoords.push_back(coordY);
                }
            }
            acceptorFile.close();
        }

        if (!donorFile.is_open()) {
            std::cerr << "No such file: " << donorConfig << "\n";
        }
        else {
                std::string line;

            while (getline(donorFile, line)) {
                if (line.empty() || line[0] == '#') {
                    continue;
                }
                else {
                    std::stringstream ss(line);

                    double coordX, coordY;
                    ss >> coordX >> coordY;

                    if (noDimension) {
                        coordX = coordX / R;
                        coordY = coordY / R;
                    }

                    donorCoords.push_back(coordX);
                    donorCoords.push_back(coordY);
                }
            }
            donorFile.close();
        }  
    }  
    
    else {
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(static_cast<long int>(now_ns));

        for (int i = 0; i < nAcceptors; ++i) {
            double randomPhi = 2.0*M_PI*randomDouble01();
            double randomR = radius*std::sqrt(randomDouble01());
            acceptorCoords.push_back(randomR*std::cos(randomPhi));
            acceptorCoords.push_back(randomR*std::sin(randomPhi));
        }

        for (int i = 0; i < nDonors; ++i) {
            double randomPhi = 2.0*M_PI*randomDouble01();
            double randomR = radius*std::sqrt(randomDouble01());
            donorCoords.push_back(randomR*std::cos(randomPhi));
            donorCoords.push_back(randomR*std::sin(randomPhi));
        }
    }
}

std::filesystem::path Configuration::getConfigFilePath(const std::string& folder, const std::string& file) {

    return std::filesystem::path(folder) / file;

}