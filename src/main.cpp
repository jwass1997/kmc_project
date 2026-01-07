#include <fstream>
#include <string>
#include <unistd.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

#include "utils.h"
#include "Configuration.h"
#include "State.h"
#include "KMCSimulator.h"

int main(int argc, char* argv[]) 
{   
    //std::cout << "Parsing arguments " << "\n";
    //argParser(argc, argv);
    int n_runs = 10;
    for (int n_r = 0; n_r < n_runs; ++n_r)
    {

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1, 1000000000);

        std::string cfg = "/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/config.txt";
        std::string acc_cfg = "/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/acceptors.txt";
        std::string don_cfg = "/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/donors.txt";
        std::string ele_cfg = "/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/electrodes.txt";
        Configuration config(cfg, acc_cfg, don_cfg, ele_cfg);
        State s0(config, dist(rng) + n_r);
        KMCSimulator kmc_alg(s0, dist(rng) + n_r);
        int n_burnin = 10000000;
        int num_bins = 100;

        int n_bin_steps = int(n_burnin / num_bins);

        std::vector<double> voltages = {0.0, 0.3, -1.2, 0.2, -0.2, 0.5, 1.5, .8};
        int o_idx = 0;
        int i_idx = 1;

        s0.updateBoundaries(voltages);
        std::string f_name = "/gpfs/bwfor/work/ws/hd_gy283-my_data/burnins_ed=0.1_a=5_steps=10e7/burnin_" + std::to_string(n_r) + ".txt";
        std::ofstream file_out(f_name);
        file_out << "t_s" << " " << "t_e" << " " << "t_mid" << " " << "I_b" << "\n";
        for (int j = 0; j < num_bins; ++j)
        {   
            std::cout << j << "\n";
            s0.resetEventCounter();
            double start_cl = s0.stateTime;;
            double current = 0.0;
            long long in_hops = 0;
            long long out_hops = 0;

            kmc_alg.simulate(s0, n_bin_steps, false, true);
            double end_cl = s0.stateTime;

            for (int k = 0; k < s0.numOfSites; ++k) 
            {
                out_hops += s0.eventCounter[(o_idx + s0.nAcceptors)*s0.numOfSites + k];
                in_hops += s0.eventCounter[s0.numOfSites*k + (o_idx + s0.nAcceptors)];
            }
            double elapsed_time = end_cl - start_cl;
            current = double(in_hops - out_hops) / elapsed_time;
            file_out << start_cl << " " << end_cl << " " << 0.5 * (start_cl + end_cl) << " " << current << "\n";
        }
        file_out.close();
    }
}