#include <random>
#include <algorithm>
#include <omp.h>
#include <boost/program_options.hpp>

#include "utils.h"

std::random_device randomDevice;
std::mt19937 rng((randomDevice()));

double sampleFromUniformDistribution(double min, double max) {

    /*
    
        Draw uniformly between min-max
    
    */

    std::uniform_real_distribution<double> uniformDist(min, max);
    double sample = uniformDist(rng);

    return sample;
}

double sampleFromNormalDistribution(double mean, double standardDeviation) {

    /*

        Returns sample from a normal distribution with mean and standardDeviation

    */

    std::normal_distribution<double> normalDistribution(mean, standardDeviation);

    double sample = normalDistribution(rng);

    return sample;
}

double calculateDistance(
    double coordinateX1, 
    double coordinateX2, 
    double coordinateY1, 
    double coordinateY2) {

    double Dx = coordinateX2 - coordinateX1;
    double Dy = coordinateY2 - coordinateY1;
    double distance = std::sqrt(Dx*Dx + Dy*Dy);
    
    return distance;
}

std::vector<double> calculateSlopes(
    std::vector<double> fX,
    std::vector<double> X) {
        
    if (fX.size() != X.size()) {
        std::cerr << "Slope can not be calculated if sizes do not match" << "\n";
    }

    std::vector<double> slopes(fX.size(), 0.0);

    double _slope = 0.0;
    for (int i = 1; i < slopes.size(); ++i) {
        _slope = (fX[i] - fX[i-1]) / (X[i] - X[i-1]);
        slopes[i] = _slope;
    }

    return slopes;
}

void createDirectoryFromStringPath(const std::string& path, const std::string& folderName) {

    /**
     * 
     * Checks if directory already exists and creates it if not. The current path is used per default.
     * 
     */

    if (folderName.empty()) {
        std::cerr << "Must specify a folder name " << "\n";
        return;
    }

    std::filesystem::path directoryPath = path.empty() ? std::filesystem::current_path() : std::filesystem::path(path);
    std::filesystem::path newFolder = std::filesystem::path(directoryPath)/folderName;

    try {
        if(std::filesystem::create_directory(newFolder)) {
            std::cout << "Folder has been created: " << newFolder << "\n";
        }
        else {
            std::cout << "Folder already exists: " << newFolder << "\n";
        }
    }
    catch(const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating new folder: " << e.what() << "\n";
    }
}