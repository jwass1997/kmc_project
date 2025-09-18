#include <iostream>
#include <vector>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <limits>

#include "Configuration.h"
#include "State.h"
#include "KMCSimulator.h"
#include "utils.h"

int main(int argc, char* argv[]) {
    argParser(argc, argv);
    /* // --- shared params ---
    const double R_nm = 150.0;                 // disk radius (nm)
    const int femRes = 100000;                   // MFEM mesh "resolution"
    const int Nr = 257;                         // lite: radial nodes
    const int Nt = 1440;                        // lite: angular nodes
    const int Nx = 401, Ny = 401;              // raster for CSV sampling
    const double xmin = -R_nm, xmax = R_nm;
    const double ymin = -R_nm, ymax = R_nm;

    const int num_electrodes = 8;
    const double width_nm = 60.0;              // arc length in nm
    const double half_width_rad = 0.5 * (width_nm / R_nm); // nm -> radians

    const std::vector<double> ang_pos_deg = {0,45,90,135,180,225,270,315};
    const std::vector<double> volts       = {-1.5,1.5,-1.5,1.5,-1.5,1.5,-1.5,1.5};

    // ======================
    // MFEM solve & export
    // ======================
    FiniteElementeCircle fem(R_nm, femRes, false);

    // define arc geometry first (0 V), then init, then set voltages, then run
    for (int k = 0; k < num_electrodes; ++k) {
        double theta0 = ang_pos_deg[k] * M_PI / 180.0;
        fem.setElectrode(0.0, theta0 - half_width_rad, theta0 + half_width_rad);
    }
    fem.initRun();
    for (int k = 0; k < num_electrodes; ++k) fem.updateElectrodeVoltage(k, volts[k]);
    fem.run();

    std::ofstream out_mfem("mfem_xy.csv");
    out_mfem << "x,y,phi\n";
    for (int iy = 0; iy < Ny; ++iy) {
        double y = ymin + (ymax - ymin) * iy / (Ny - 1);
        for (int ix = 0; ix < Nx; ++ix) {
            double x = xmin + (xmax - xmin) * ix / (Nx - 1);
            if (std::hypot(x,y) <= R_nm + 1e-9) {
                out_mfem << x << "," << y << "," << fem.getPotential(x,y) << "\n";
            }
        }
    }
    std::cout << "Wrote mfem_xy.csv\n";

    // ======================
    // Lite solve & export
    // ======================
    LiteLaplaceCircle solver(R_nm, Nr, Nt, true, 1.85);

    for (int k = 0; k < num_electrodes; ++k) {
        double center_rad = ang_pos_deg[k] * M_PI / 180.0;
        solver.setElectrode(volts[k], center_rad - half_width_rad, center_rad + half_width_rad);
    }

    int iters = solver.run(50000, 1e-10);
    std::cout << "Lite solver converged in " << iters << " iterations\n";

    std::ofstream out_lite("lite_xy.csv");
    out_lite << "x,y,phi\n";
    for (int iy = 0; iy < Ny; ++iy) {
        double y = ymin + (ymax - ymin) * iy / (Ny - 1);
        for (int ix = 0; ix < Nx; ++ix) {
            double x = xmin + (xmax - xmin) * ix / (Nx - 1);
            if (std::hypot(x,y) <= R_nm + 1e-9) {
                out_lite << x << "," << y << "," << solver.getPotential(x,y) << "\n";
            }
        }
    }
    std::cout << "Wrote lite_xy.csv\n";

    return 0; */
}