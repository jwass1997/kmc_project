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
    argParser(argc, argv);
    /*int n_runs = 10;
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
        int n_burnin = 1000000;
        int num_bins = 100;

        int n_bin_steps = int(n_burnin / num_bins);

        std::vector<double> voltages = {0.0, 0.3, -1.2, 0.2, -0.2, 0.5, .5, .8};
        int o_idx = 0;
        int i_idx = 1;

        s0.updateBoundaries(voltages);
        std::string f_name = "/gpfs/bwfor/work/ws/hd_gy283-my_data/burnins_ed=0.01_a=20_steps=10e6/burnin_" + std::to_string(n_r) + ".txt";
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
    }*/
}

/*#include <fstream>
#include <string>
#include <unistd.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <random>
#include <limits>
#include <iomanip>
#include <filesystem>

#include "utils.h"
#include "Configuration.h"
#include "State.h"
#include "KMCSimulator.h"

int main(int argc, char* argv[])
{
    argParser(argc, argv);*/

    /* int n_runs = 10;

    for (int n_r = 0; n_r < n_runs; ++n_r)
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1, 1000000000);

        std::string cfg     = "/home/hd/hd_hd/hd_gy283/kmc_project/configs/std_configs/config.txt";
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
        int o_idx = 0; // your electrode/site index convention
        int i_idx = 1; // unused here, kept from your code

        s0.updateBoundaries(voltages);

        // Output per-bin current
        std::string out_dir = "/gpfs/bwfor/work/ws/hd_gy283-my_data/burnins_ed=0.01_a=20_steps=10e7/";
        std::error_code ec;
        std::filesystem::create_directories(out_dir, ec); // creates parents too
        if (ec) {
            std::cerr << "Failed to create directory '" << out_dir << "': " << ec.message() << "\n";
        }
        std::string f_name = out_dir + "burnin_" + std::to_string(n_r) + ".txt";
        std::ofstream file_out(f_name);
        file_out << "t_s" << " " << "t_e" << " " << "t_mid" << " " << "I_b" << "\n";

        // Store bin durations and currents for blocking analysis
        std::vector<double> dt_bins(num_bins, 0.0);
        std::vector<double> I_bins(num_bins, 0.0);

        for (int j = 0; j < num_bins; ++j)
        {
            std::cout << j << "\n";
            s0.resetEventCounter();
            double start_cl = s0.stateTime;
            long long in_hops = 0;
            long long out_hops = 0;

            kmc_alg.simulate(s0, n_bin_steps, false, true);
            double end_cl = s0.stateTime;

            for (int k = 0; k < s0.numOfSites; ++k)
            {
                out_hops += s0.eventCounter[(o_idx + s0.nAcceptors) * s0.numOfSites + k];
                in_hops  += s0.eventCounter[s0.numOfSites * k + (o_idx + s0.nAcceptors)];
            }

            double elapsed_time = end_cl - start_cl;
            double current = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(elapsed_time) && elapsed_time > 0.0) {
                current = double(in_hops - out_hops) / elapsed_time;
            }

            double t_mid = 0.5 * (start_cl + end_cl);
            file_out << start_cl << " " << end_cl << " " << t_mid << " " << std::setprecision(12) << current << "\n";

            dt_bins[j] = elapsed_time;
            I_bins[j]  = current;
        }
        file_out.close();

        // ---- Blocking analysis: compute SE(L) for many block sizes L ----
        // We use TIME-WEIGHTED block means because bin durations vary:
        //   I_block = sum(I_j * dt_j) / sum(dt_j)
        //
        // Then: SE(L) = std(block_means, ddof=1) / sqrt(B_effective)
        // where B_effective is the number of valid blocks (with sum(dt)>0).

        std::string se_name = out_dir + "blocking_se_" + std::to_string(n_r) + ".txt";
        std::ofstream se_out(se_name);
        se_out << "L_bins" << " " << "B_blocks" << " " << "SE_timeweighted" << "\n";

        auto finite_pos = [](double x) { return std::isfinite(x) && x > 0.0; };

        int min_B = 10; // require at least this many blocks to trust the estimate
        int max_L = std::max(1, num_bins / min_B);

        for (int L = 1; L <= max_L; L *= 2) // powers of 2: 1,2,4,8,...
        {
            int B = num_bins / L;
            if (B < 2) continue;

            std::vector<double> block_means;
            block_means.reserve(B);

            for (int b = 0; b < B; ++b)
            {
                double sum_dt = 0.0;
                double sum_Idt = 0.0;

                int start = b * L;
                int end   = start + L;

                for (int j = start; j < end; ++j)
                {
                    if (!finite_pos(dt_bins[j]) || !std::isfinite(I_bins[j])) continue;
                    sum_dt  += dt_bins[j];
                    sum_Idt += I_bins[j] * dt_bins[j];
                }

                if (sum_dt > 0.0) {
                    block_means.push_back(sum_Idt / sum_dt);
                }
            }

            int Beff = static_cast<int>(block_means.size());
            if (Beff < 2) {
                se_out << L << " " << Beff << " " << std::numeric_limits<double>::quiet_NaN() << "\n";
                continue;
            }

            // sample mean
            double mean = 0.0;
            for (double x : block_means) mean += x;
            mean /= static_cast<double>(Beff);

            // sample variance (ddof=1)
            double var = 0.0;
            for (double x : block_means) {
                double d = x - mean;
                var += d * d;
            }
            var /= static_cast<double>(Beff - 1);

            double se = std::sqrt(var) / std::sqrt(static_cast<double>(Beff));
            se_out << L << " " << Beff << " " << std::setprecision(12) << se << "\n";
        }

        se_out.close();
    }

    return 0; */
//}
