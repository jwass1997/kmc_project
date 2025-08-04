#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>
#include <fstream>

#include "Configuration.h"
#include "State.h"
#include "KMCSimulator.h"
#include "FEMmethods.h"
#include "utils.h"
#include "Random.h"

int main(int argc, char* argv[]) {
    //argParser(argc, argv);

    /* std::string cfg = "/gpfs/bwfor/home/hd/hd_hd/hd_gy283/kmc_project/configs/config.txt";
    std::string acc = "/gpfs/bwfor/home/hd/hd_hd/hd_gy283/kmc_project/configs/acceptors.txt";
    std::string don = "/gpfs/bwfor/home/hd/hd_hd/hd_gy283/kmc_project/configs/donors.txt";
    std::string ele = "/gpfs/bwfor/home/hd/hd_hd/hd_gy283/kmc_project/configs/electrodes.txt";

    Configuration config(cfg, acc, don, ele, false);

    State state(config);

    std::vector<double> voltages = {1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    state.updateBoundaries(voltages);
    KMCSimulator kmc(state);
    kmc.simulate(state, 10000, false, true);
    for (int i = 0; i < state.numOfSites*state.numOfSites; ++i) {
        std::cout << state.eventCounter[i];
    } */
    //std::cout << state.currentPotential[200] << "\n";
    //State state2(state);
    //std::cout << state2.currentPotential[200] << "\n";
    //KMCSimulator kmc(state);
    
    /* std::ofstream ofs("phi_0.txt");
    if (!ofs) {
        std::cerr << "Cannot open output file\n";
        return 1;
    }

    const int radialSteps = 100;
    const int angularSteps = 360;
    for (int ir = 1; ir <= radialSteps; ++ir) {
        double r = config.radius * ir / radialSteps;
        for (int jt = 0; jt < angularSteps; ++jt) {
            double theta = 2 * M_PI * jt / angularSteps;
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            double phi = fem.getPotential(x, y);
            ofs << x << " " << y << " " << phi << "\n";
        }
        ofs << "\n";
    }

    std::vector<double> newBoundaries = {10.5, -9.2, 9.9, 8.2, -9.0, 8.5, 9.4, -9.1};
    state.updateBoundaries(newBoundaries, fem);

    std::ofstream ofs_1("phi_1.txt");
    if (!ofs_1) {
        std::cerr << "Cannot open output file\n";
        return 1;
    }
    for (int ir = 1; ir <= radialSteps; ++ir) {
        double r = config.radius * ir / radialSteps;
        for (int jt = 0; jt < angularSteps; ++jt) {
            double theta = 2 * M_PI * jt / angularSteps;
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            double phi = fem.getPotential(x, y);
            ofs_1 << x << " " << y << " " << phi << "\n";
        }
        ofs_1 << "\n";
    }

    return 0; */
}