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
    argParser(argc, argv);
}

/* #include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include "FEMmethods.h" // your original header

// Safe fallback for PI if not provided
#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

// Sample the circular domain uniformly in polar coords and dump potentials to CSV.
// samplesR: radial subdivisions (including center), samplesTheta: angular subdivisions.
void dumpCircleSamplingCSV(FiniteElementeCircle &fe,
                           double radius,
                           int samplesR,
                           int samplesTheta,
                           const std::string &fname) {
    std::ofstream ofs(fname);
    if (!ofs) {
        std::cerr << "Failed to open " << fname << " for writing\n";
        return;
    }
    ofs << "x,y,potential\n";
    for (int ir = 0; ir < samplesR; ++ir) {
        double r = radius * ir / double(samplesR - 1); // 0..radius
        for (int it = 0; it < samplesTheta; ++it) {
            double theta = 2.0 * M_PI * it / double(samplesTheta);
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            double pot = NAN;
            try {
                pot = fe.getPotential(x, y);
            } catch (const std::exception &e) {
                // out-of-range should not happen if r <= radius; leave as NaN
            }
            ofs << x << "," << y << "," << pot << "\n";
        }
    }
}

int main() {
    try {
        // domain / mesh parameters
        double radius = 1.0;
        int maxElements = 2000; // adjust for resolution
        bool saveSolution = true;

        // sampling resolution for output CSV
        const int samplesR = 80;
        const int samplesTheta = 120;

        // Construct the circular finite element problem
        FiniteElementeCircle fe(radius, maxElements, saveSolution);
        fe.initRun();

        // First solve: electrode from angle 0 to pi/4 at voltage 1.0
        double begin_angle = 0.0;
        double end_angle = M_PI / 4.0;
        fe.setElectrode(1.0, begin_angle, end_angle);
        fe.run(); // solution #0
        dumpCircleSamplingCSV(fe, radius, samplesR, samplesTheta, "circle_electrode_1.0.csv");

        // Update electrode voltage to 2.0 and resolve
        fe.updateElectrodeVoltage(0, 2.0); // assuming this is the first electrode
        fe.run(); // solution #1
        dumpCircleSamplingCSV(fe, radius, samplesR, samplesTheta, "circle_electrode_2.0.csv");

        std::cout << "Done. Outputs:\n"
                  << "  circle_electrode_1.0.csv\n"
                  << "  circle_electrode_2.0.csv\n"
                  << "  (and GF files from original run: laplace_solution*.gf)\n";
    } catch (const std::exception &e) {
        std::cerr << "Exception during run: " << e.what() << std::endl;
        return 1;
    }
    return 0;
} */


