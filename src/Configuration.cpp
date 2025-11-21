#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <queue>

#include "Configuration.h"

bool subset_connected(
    std::vector<double> &subset_coords, 
    std::vector<double> &coords, 
    double min_dist, 
    double max_dist)
    {
        if (min_dist > max_dist)
        {
            throw std::runtime_error("Minimal hopping distance is greater than maximal hopping distance.");
        }

        const int N_main = static_cast<int>(coords.size() / 2);
        const int N_subset = static_cast<int>(subset_coords.size() / 2);

        if (N_subset <= 1) 
        {
            return true;
        }

        std::vector<double> nodes;
        nodes.reserve(coords.size() + subset_coords.size());
        nodes.insert(nodes.end(), subset_coords.begin(), subset_coords.end());
        nodes.insert(nodes.end(), coords.begin(), coords.end());

        const int total_points = static_cast<int>(nodes.size() / 2);
        const int subset_start = 0;
        const int subset_end = N_subset;

        const double mind2 = min_dist * min_dist;
        const double maxd2 = max_dist * max_dist;

        std::vector<std::vector<int>> adj(total_points);

        for (int i = 0; i < total_points; ++i)
        {   
            double x_i = nodes[2*i];
            double y_i = nodes[2*i+1];

            for (int j = i+1; j < total_points; ++j)
            {
                double dx = x_i - nodes[2*j];
                double dy = y_i - nodes[2*j+1];

                double dist2 = dx*dx + dy*dy;
                
                if (dist2 >= mind2 && dist2 <= maxd2)
                {
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                }
            }

            if (i < N_subset && adj[i].empty())
            {
                return false;
            }
        }

        std::vector<bool> visited(total_points, false);
        std::queue<int> q;

        int start = subset_start;
        visited[start] = true;
        q.push(start);

        while (!q.empty())
        {
            int u = q.front();
            q.pop();

            for (int v : adj[u])
            {
                if (!visited[v])
                {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        for (int idx = subset_start; idx <  subset_end; ++idx)
        {
            if (!visited[idx])
            {
                return false;
            }
        }

        return true;
    }

std::vector<double> gm_rej_sampling(
    int num_samples, 
    double radius,
    double min_distance,
    double max_distance,
    GaussianMixture2D& gaussian_mixture,
    std::vector<double>& subset_nodes,
    const int max_attempts) {

        if (max_distance < min_distance)
        {
            throw std::runtime_error("Minimal distance > maximal distance.");
        }

        double radius2 = radius*radius;

        std::vector<double> node_samples;
        for  (int attempt = 0; attempt < max_attempts; ++attempt) 
        {
            node_samples.clear();
            while ((int)node_samples.size() < 2 * num_samples)
            {
                std::array<double, 2> pnt = gaussian_mixture.sample();

                if (pnt[0]*pnt[0] + pnt[1]*pnt[1] < radius2)
                {
                    node_samples.push_back(pnt[0]);
                    node_samples.push_back(pnt[1]);
                }
                else 
                {
                    continue;            
                }
            }
            
            bool connected_nodes = subset_connected(subset_nodes, node_samples, min_distance, max_distance);

            if (connected_nodes) 
            {
                return node_samples;            
            }
        }
        
        return {};
}

Configuration::Configuration(ConfigurationParams& params, uint64_t seed)
    : _rng(seed)
{

    if (params.nElectrodes != params.electrodeData.size()) 
    {
        throw std::invalid_argument("[Configuration]: electrodeData.size() is not equal to nElectrodes");
    }

    nAcceptors = params.nAcceptors;
    nDonors = params.nDonors;
    nElectrodes = params.nElectrodes;

    numOfSites = nAcceptors + nElectrodes;

    radius = params.radius;

    nu0 = params.nu0;
    a = params.a;
    T = params.T;

    kbT = kb*T;

    energyDisorder = params.energyDisorder * e / (kb*T);

    electrodeWidth = params.electrodeWidth;
    minHopDistance = params.minHopDistance;
    maxHopDistance = params.maxHopDistance;

    noDimension = params.noDimension;

    Nr = params.Nr;
    Nt = params.Nt;

    A0 = (e*e) / (4.0*kb*T*PI*eps0*epsr*1e-9);
    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

    if (noDimension) {
        A0 = A0 / R;
        a = a / R;
        radius = radius / R;
        electrodeWidth = electrodeWidth / R;
        minHopDistance = minHopDistance / R;
        maxHopDistance = maxHopDistance / R;        
    }

    distType = params.distType;

    electrodeData = params.electrodeData;
    electrodeCoords.reserve(2*nElectrodes);
    for (const auto& el : electrodeData) 
    {
        double phi = (2.0 * M_PI * el.angularPosition) / 360.0;
        double x = radius * std::cos(phi);
        double y = radius * std::sin(phi);
        electrodeCoords.push_back(x);
        electrodeCoords.push_back(y);
    }

    if (distType == "uniform") 
    {
        for (int i = 0; i < nAcceptors; ++i) 
        {
            double randomPhi = 2.0*M_PI*uniform01();
            double randomR = radius*std::sqrt(uniform01());
            acceptorCoords.push_back(randomR*std::cos(randomPhi));
            acceptorCoords.push_back(randomR*std::sin(randomPhi));
        }

        for (int i = 0; i < nDonors; ++i) 
        {
            double randomPhi = 2.0*M_PI*uniform01();
            double randomR = radius*std::sqrt(uniform01());
            donorCoords.push_back(randomR*std::cos(randomPhi));
            donorCoords.push_back(randomR*std::sin(randomPhi));
        }
    }

    if (distType == "gaussian_mixture")
    {
        const int max_attempts = 100;

        double std_min = radius / 10.0;
        double std_max = radius / 2.0;

        double radius_min = 0.0;
        double radius_max = 0.7 * radius;

        //double std_min = std::max(minHopDistance / 3.0, 0.02*radius);
        //double std_max = std::min(maxHopDistance / 2.0, 0.5*radius);

        int num_components = params.n_comps;

        while (acceptorCoords.empty())
        {
            std::vector<Gaussian2D> mixture_components;
            std::vector<double> component_weights;

            for (int i = 0; i < num_components; ++i)
            {
                double comp_weight = uniform01();
                component_weights.push_back(comp_weight);

                double randomPhi = 2.0*M_PI*uniform01();
                double randomR = (radius_max - radius_min) * std::sqrt(uniform01()) + radius_min;

                double mx = randomR*std::cos(randomPhi);
                double my = randomR*std::sin(randomPhi);

                double log_std_min = std::log(std_min); 
                double log_std_max = std::log(std_max);
                
                double log_sx = (log_std_max - log_std_min) * uniform01() + log_std_min;
                double log_sy = (log_std_max - log_std_min) * uniform01() + log_std_min;

                double sx = std::exp(log_sx);
                double sy = std::exp(log_sy);

                //double sx = (std_max - std_min) * uniform01() + std_min;
                //double sy = (std_max - std_min) * uniform01() + std_min;

                double varxx = sx * sx;
                double varyy = sy * sy;
                double varxy = 0.0;

                mixture_components.push_back(Gaussian2D(mx, my, 0.0, varxx, varyy));

                params.var11.push_back(varxx);
                params.var22.push_back(varyy);
                params.var12.push_back(varxy);

                params.m1.push_back(mx);
                params.m2.push_back(my);
            }

            GaussianMixture2D gaussian_mixture(mixture_components, component_weights);

            params.normalized_weights = gaussian_mixture.normalized_weights;

            acceptorCoords = gm_rej_sampling(nAcceptors, radius, minHopDistance, maxHopDistance, gaussian_mixture, electrodeCoords, max_attempts);
        }

        for (int k = 0; k < nDonors; ++k) 
        {
            double randomPhi = 2.0*M_PI*uniform01();
            double randomR = radius*std::sqrt(uniform01());
            donorCoords.push_back(randomR*std::cos(randomPhi));
            donorCoords.push_back(randomR*std::sin(randomPhi));
        }
    }
}

Configuration::Configuration(
    const std::string& cfg, 
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg
) {

    auto config = cfg;
    auto acceptorConfig = acceptorCfg;
    auto donorConfig = donorCfg;
    auto electrodeConfig = electrodeCfg;
    /* std::cout
    << "  CWD:    " << std::filesystem::current_path() << "\n"
    << "  Target: " << std::filesystem::absolute(config) << "\n"
    << "  Exists: " << std::boolalpha << std::filesystem::exists(config) << "\n"; */
    std::ifstream configFile(config);
    std::ifstream acceptorFile(acceptorConfig);
    std::ifstream donorFile(donorConfig);
    std::ifstream electrodeFile(electrodeConfig);
    //std::cout << config << "\n";
    if (!configFile.is_open()) {
        std::cerr << "[FATAL] Could not open config file: " << config << "\n";
        std::exit(1);
    }
    if (!acceptorFile.is_open()) {
        std::cerr << "[FATAL] Could not open acceptor config: " << acceptorConfig << "\n";
        std::exit(1);
    }
    if (!donorFile.is_open()) {
        std::cerr << "[FATAL] Could not open donor config: " << donorConfig << "\n";
        std::exit(1);
    }
    if (!electrodeFile.is_open()) {
        std::cerr << "[FATAL] Could not open electrode config: " << electrodeConfig << "\n";
        std::exit(1);
    }

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
                else if (key == "Nr") {
                    Nr = std::stoi(value);
                }
                else if (key == "Nt") {
                    Nt = std::stoi(value);
                }
            }
        }
        configFile.close();
    }

    if (nAcceptors <= 0) 
    {
    std::cerr << "[FATAL] Parsed invalid nAcceptors: " << nAcceptors << "\n";
    std::exit(1);
    }
    if (nElectrodes <= 0) 
    {
        std::cerr << "[FATAL] Parsed invalid nElectrodes: " << nElectrodes << "\n";
        std::exit(1);
    }

    numOfSites = nAcceptors + nElectrodes;
    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));
    kbT = kb*T;
    A0 = (e*e) / (4.0*kb*T*PI*eps0*epsr*1e-9);

    if (noDimension) 
    {
        radius = radius / R;
        a = a / R;
        A0 = A0 / R;
        electrodeWidth = electrodeWidth / R;
        minHopDistance = minHopDistance / R;
        maxHopDistance = maxHopDistance / R;
    }

    /* Electrodes */
    if (!electrodeFile.is_open()) 
    {
            std::cerr << "No such file: " << electrodeConfig<< "\n";
    }
    
    else 
    {
        std::string line;
        
        while (getline(electrodeFile, line)) 
        {
            if (line.empty() || line[0] == '#') 
            {
                continue;
            }
            else 
            {
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
    /* Coordinates */
    if (!acceptorFile.is_open()) 
    {
        std::cerr << "No such file: " << acceptorConfig << "\n";
    }
    else 
    {
        std::string line;

        while (getline(acceptorFile, line)) 
        {
            if (line.empty() || line[0] == '#') 
            {
                continue;
            }
            else 
            {
                std::stringstream ss(line);

                double coordX, coordY;
                ss >> coordX >> coordY;
                
                if (noDimension) 
                {
                    coordX = coordX / R;
                    coordY = coordY / R;
                }

                acceptorCoords.push_back(coordX);
                acceptorCoords.push_back(coordY);
            }
        }
        acceptorFile.close();
    }

    if (!donorFile.is_open()) 
    {
        std::cerr << "No such file: " << donorConfig << "\n";
    }
    else 
    {
            std::string line;

        while (getline(donorFile, line)) 
        {
            if (line.empty() || line[0] == '#') 
            {
                continue;
            }
            else 
            {
                std::stringstream ss(line);

                double coordX, coordY;
                ss >> coordX >> coordY;

                if (noDimension) 
                {
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

double Configuration::uniform01() {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    return U(_rng);
}