#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>

#include "Configuration.h"
#include "State.h"
#include "KMCSimulator.h"
#include "FEMmethods.h"
#include "utils.h"
#include "Random.h"

int main() {

    setRandomSeed(42);

    std::string configsPath = "configs";
    std::string dataPath = "data";

    recordDevice("1", 1e4, 1e5, configsPath, dataPath);
}