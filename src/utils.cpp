#include <random>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <boost/program_options.hpp>

#include "utils.h"
#include "Random.h"
#include "State.h"
#include "FEMmethods.h"
#include "Configuration.h"
#include "KMCSimulator.h"

double calculateDistance(
    double coordinateX1, 
    double coordinateX2, 
    double coordinateY1, 
    double coordinateY2
) {

    double Dx = coordinateX2 - coordinateX1;
    double Dy = coordinateY2 - coordinateY1;
    double distance = std::sqrt(Dx*Dx + Dy*Dy);
    
    return distance;
}

std::vector<std::vector<double>> scaledLHC(
    std::size_t N,
    std::size_t D,
    const std::vector<double>& mins,
    const std::vector<double>& maxs,
    unsigned int seed
) {
    if (mins.size() != D || maxs.size() != D)
        throw std::invalid_argument("mins/maxs must have length D");

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // Step 1: generate in [0,1]^D
    std::vector<std::vector<double>> samples(N, std::vector<double>(D));
    for (std::size_t d = 0; d < D; ++d) {
        std::vector<std::size_t> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng);

        for (std::size_t i = 0; i < N; ++i) {
            double u = uni(rng);
            samples[i][d] = (idx[i] + u) / static_cast<double>(N);
        }
    }

    // Step 2: scale to [mins[d], maxs[d]]
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t d = 0; d < D; ++d) {
            samples[i][d] = mins[d] + samples[i][d] * (maxs[d] - mins[d]);
        }
    }

    return samples;
}

void singleRun(
    const std::string& ID, 
    int equilibriumSteps, 
    int numOfSteps, 
    const std::string& cfg, 
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolderPath
) {

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    int seed0 = 1234567890;
    auto now = std::chrono::high_resolution_clock::now();
    auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    setRandomSeed(seed0 + static_cast<long int>(now_ns));

    Configuration config(cfg, acceptorCfg, donorCfg, electrodeCfg, false);

    State state(config);
    KMCSimulator simulator(state);

    int nAcceptors = state.nAcceptors;
    int nElectrodes = state.nElectrodes;
    int nDonors = state.nDonors;

    std::vector<double> flattenedAcceptorCoordinates;
    std::vector<double> flattenedDonorCoordinates;
    std::vector<double> flattenedElectrodeCoordinates;    
    std::vector<int> flattenedEventCounts;

    std::vector<size_t> shapeFlattenedAcceptorCoordinates = {static_cast<size_t>(nAcceptors), 2};
    std::vector<size_t> shapeFlattenedDonorCoordinates = {static_cast<size_t>(nDonors), 2};
    std::vector<size_t> shapeFlattenedElectrodeCoordinates = {static_cast<size_t>(nElectrodes), 2};
    std::vector<size_t> shapeFlattenedEventCounts = {static_cast<size_t>(nAcceptors+nElectrodes), static_cast<size_t>(nAcceptors+nElectrodes)};

    std::string deviceName = saveFolderPath + "/device_" + ID + ".npz";
    cnpy::npz_save(deviceName, "ID", &ID, {1}, "w"); 

    simulator.simulate(state, equilibriumSteps, true, false);
    simulator.simulate(state, numOfSteps, false, true);

    for(int i = 0; i < nAcceptors; ++i) {
        flattenedAcceptorCoordinates.push_back(state.acceptorCoordinates[i*2]);
        flattenedAcceptorCoordinates.push_back(state.acceptorCoordinates[i*2 + 1]);
    }
    for(int i = 0; i < nDonors; ++i) {
        flattenedDonorCoordinates.push_back(state.donorCoordinates[i*2]);
        flattenedDonorCoordinates.push_back(state.donorCoordinates[i*2 + 1]);
    }
    for(int i = 0; i < nElectrodes; ++i) {
        flattenedElectrodeCoordinates.push_back(state.electrodeCoordinates[i*2]);
        flattenedElectrodeCoordinates.push_back(state.electrodeCoordinates[i*2 + 1]);
    }

    for(int j = 0; j < nAcceptors+nElectrodes; ++j) {
        for(int i = 0; i <nAcceptors+nElectrodes; ++i) {
            flattenedEventCounts.push_back(state.eventCounter[j*state.numOfSites + i]);
        }
    }

    double total_time = state.stateTime;

    cnpy::npz_save(deviceName, "acceptor_coordinates", flattenedAcceptorCoordinates.data(), shapeFlattenedAcceptorCoordinates, "a");
    cnpy::npz_save(deviceName, "donor_coordinates", flattenedDonorCoordinates.data(), shapeFlattenedDonorCoordinates, "a");
    cnpy::npz_save(deviceName, "electrode_coordinates", flattenedElectrodeCoordinates.data(), shapeFlattenedElectrodeCoordinates, "a");
    cnpy::npz_save(deviceName, "event_counts", flattenedEventCounts.data(), shapeFlattenedEventCounts, "a");
    cnpy::npz_save(deviceName, "device_time", &total_time, {1}, "a");
}

double singleIVPoint(
    State& initState,
    int outputIdx,
    int numOfTasks,
    int simSteps,
    int femRes,
    std::vector<double> voltages
) {
    State equilState(initState); 
    KMCSimulator kmc(equilState);

    equilState.resetEventCounter();

    double averagedCurrent = 0.0;
    double totalTime = 0.0;
    int intervalSteps = simSteps / numOfTasks;
    int netEvents = 0;

    int intervalCount = 0;
    while (intervalCount < numOfTasks) {

        double startClock = equilState.stateTime;
        kmc.simulate(equilState, intervalSteps, false, true);
        double endClock = equilState.stateTime; 

        double elapsedTime = endClock - startClock;
        int inEvents = 0;
        int outEvents = 0;
        for (int i = 0; i < equilState.numOfSites; ++i) {
            outEvents += equilState.eventCounter[(outputIdx + equilState.nAcceptors)*equilState.numOfSites + i];
            inEvents += equilState.eventCounter[equilState.numOfSites*i + (outputIdx + equilState.nAcceptors)];
        }
        totalTime += elapsedTime;
        netEvents += inEvents-outEvents;

        equilState.resetEventCounter();

        intervalCount++;
    }

    averagedCurrent = static_cast<double>(netEvents) / totalTime;

    return averagedCurrent;    
}

void singleIVCurve(
    int numOfPoints,
    int inputIdx,
    int outputIdx,
    double minVoltage,
    double maxVoltage,
    int eqSteps,
    int simSteps,
    int numOfTasks,
    int seed,
    std::vector<double> controlVoltages,
    const std::string& cfg,
    const std::string& accCfg,
    const std::string& donCfg,
    const std::string& eleCfg,
    const std::string& saveFolder,
    const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("singleIVCurve(): Save folder not found");
    }

    std::string file = saveFolder + "/" + fileName + ".npz";

    std::vector<double> currentData(numOfPoints, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(numOfPoints)};

    std::vector<double> controlData(controlVoltages.size(), 0.0);
    std::vector<size_t> controlDataShape = {controlVoltages.size()};

    for (int cv = 0; cv < controlVoltages.size(); ++cv) {
        controlData[cv] = controlVoltages[cv];
    }
    controlData[inputIdx] = -999.999;
    controlData[outputIdx] = -999.999;

    int femRes = 100000;

    double range = maxVoltage - minVoltage;
    double vStep = range / static_cast<double>(numOfPoints);

    Configuration config(
        cfg,
        accCfg,
        donCfg,
        eleCfg,
        false
    );
    State state(config);
    KMCSimulator kmc(state);
    kmc.simulate(state, eqSteps, false, false);
    state.resetEventCounter();
    state.stateTime = 0.0;

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int v = 0; v < numOfPoints; ++v) {
            
            std::vector<double> voltages(controlVoltages.size(), 0.0);
            voltages[outputIdx] = 0.0;
            voltages[inputIdx] = minVoltage + v*vStep;
            singleIVPoint(
                state,
                outputIdx,
                numOfTasks,
                simSteps,
                femRes,
                voltages        
            );
        }
    }
}

int argParser(int argc, char* argv[]) {

    std::cout << "Hello argParser!" << "\n";

    boost::program_options::options_description globalOptions(" ");

    globalOptions.add_options()
        ("help, h", "help message")
        ("command", boost::program_options::value<std::string>(), "command to run")
    ;

    boost::program_options::positional_options_description position;
    position.add("command", 1);

    auto parseCommand = boost::program_options::command_line_parser(argc, argv).options(globalOptions).positional(position).allow_unregistered().run();

    boost::program_options::variables_map commandVM;
    boost::program_options::store(parseCommand, commandVM);

    if (commandVM.count("help") || !commandVM.count("command")) {
        std::cout << globalOptions << "\n"
                  << "Allowed commands:\n"
                  << " singleRun --configs <string> --save_path <string> --equilibriumSteps <int> --simulationSteps <int> --deviceName <string>\n";
        return 0;
    }

    std::string firstCommand = commandVM["command"].as<std::string>();
    std::vector<std::string> remainingCommand = boost::program_options::collect_unrecognized(
        parseCommand.options,
        boost::program_options::include_positional
    );

    remainingCommand.erase(remainingCommand.begin());

    if (firstCommand == "singleRun") {

        boost::program_options::options_description options("Single run options");
        options.add_options()
            ("cfg", boost::program_options::value<std::string>()->required())
            ("acceptorCfg", boost::program_options::value<std::string>()->required())
            ("donorCfg", boost::program_options::value<std::string>()->required())
            ("electrodeCfg", boost::program_options::value<std::string>()->required())
            ("save_path", boost::program_options::value<std::string>()->required())
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("deviceName", boost::program_options::value<std::string>()->required())
        ;

        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        singleRun(
            vm["deviceName"].as<std::string>(),
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["save_path"].as<std::string>()
        );

        return 1;
    }

    return 1;
}