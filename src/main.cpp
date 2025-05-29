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

int main(int argc, char* argv[]) {

    std::string configsPath = "../configs";
    std::string dataPath = "../data";
    /* singleStateBatch(
        1000,
        0,
        -1.5,
        1.5,
        1e4, 
        1e5, 
        100, 
        configsPath,
        dataPath,
        "testing"
    );
    singleRun("1", 1e4, 1e5, configsPath, dataPath); */

    argParser(argc, argv);
}