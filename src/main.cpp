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

inline double wrap_angle(double x) {
    x = std::fmod(x, 2.0 * PI);
    if (x < 0) x += 2.0 * PI;
    return x;
}

double von_mises_sample(double mu, double kappa, std::mt19937& gen) {
    std::uniform_real_distribution<double> U(0.0, 1.0);

    if (kappa < 1e-12) {
        return wrap_angle(mu + U(gen) * 2.0 * PI);
    }

    double a = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    double b = (a - std::sqrt(2.0 * a)) / (2.0 * kappa);
    double r = (1.0 + b * b) / (2.0 * b);

    while (true) {
        double u1 = U(gen);
        double u2 = U(gen);

        double z  = std::cos(PI * u1);
        double f  = (1.0 + r * z) / (r + z);
        double c  = kappa * (r - f);

        if (u2 < c * (2.0 - c) || u2 <= c * std::exp(1.0 - c)) {
            double u3 = U(gen);
            double theta = (u3 < 0.5 ? -1.0 : 1.0) * std::acos(f) + mu;
            return wrap_angle(theta);
        }
    }
}

double beta_sample(double alpha, double beta, std::mt19937& gen) {
    if (alpha <= 0.0 || beta <= 0.0) throw std::invalid_argument("alpha,beta > 0");
    std::gamma_distribution<double> G1(alpha, 1.0);
    std::gamma_distribution<double> G2(beta, 1.0);
    double x = G1(gen);
    double y = G2(gen);
    return x / (x + y); // in (0,1)
}

int main(int argc, char* argv[]) {

    /* std::mt19937 rng(412412312);

    std::ofstream file;
    int n_a = 200;
    double R = 150.0;
    double mu = 0.0;
    double kappa = 0.2;
    double alpha = 4.5;
    double beta = 4.5;
    file.open("/home/hd/hd_hd/hd_gy283/kmc_project/cpp_vMB.txt");
    int counter = 0;
    while (counter < n_a) {
        double u = beta_sample(alpha, beta, rng);
        double r = R*std::sqrt(u);        
        double angle = von_mises_sample(mu, kappa, rng);

        double x = r*std::cos(angle);
        double y = r*std::sin(angle);

        if (std::sqrt(x*x + y*y) < R) {
            file << x << " " << y << "\n";
            ++counter;
        }            
    }
    file.close(); */

    argParser(argc, argv);

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