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
    int eqSteps, 
    int simSteps, 
    std::vector<double> voltages,
    const std::string& cfg, 
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    int seed,
    const std::string& saveFolderPath
) {

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("[singleRun]: No save folder specified !");
    }

    setRandomSeed(seed);

    Configuration config(
        cfg, 
        acceptorCfg, 
        donorCfg, 
        electrodeCfg
    );

    State state(config);
    KMCSimulator kmc(state);
    kmc.simulate(state, eqSteps, false, false);
    state.resetEventCounter();
    state.stateTime = 0.0;

    state.updateBoundaries(voltages);
    kmc.simulate(state, simSteps, false, true);

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

    std::string deviceName = saveFolderPath + "/" + ID + ".npz";
    cnpy::npz_save(deviceName, "ID", &ID, {1}, "w"); 

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
    std::vector<double> voltages
) {
    State equilState(initState);
    equilState.updateBoundaries(voltages);

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
        throw std::invalid_argument("[singleIVCurve]: Save folder not found");
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

    double range = maxVoltage - minVoltage;
    double vStep = range / static_cast<double>(numOfPoints);

    Configuration config(
        cfg,
        accCfg,
        donCfg,
        eleCfg
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
            for (int k = 0; k < voltages.size(); ++k) {
                voltages[k] = controlVoltages[k];
                //std::cout << "VOLTAGE " << k << " " << voltages[k] << "\n";
            }
            voltages[outputIdx] = 0.0;
            voltages[inputIdx] = minVoltage + v*vStep;
            double currentOutput = singleIVPoint(
                state, 
                outputIdx, 
                numOfTasks, 
                simSteps, 
                voltages
            );

            currentData[v] = currentOutput;
        }
    }

    cnpy::npz_save(file, "current", currentData.data(), currentDataShape, "w");
    cnpy::npz_save(file, "control", controlData.data(), controlDataShape, "a");
    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "a");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
}

void batchOfIVPoints(
    int batchSize,
    double minVoltage,
    double maxVoltage,
    int inputIdx,
    int outputIdx,
    int eqSteps,
    int simSteps,
    int numOfTasks,
    int LHCSeed,
    int threadBaseSeed,
    const std::string& cfg,
    const std::string& accCfg,
    const std::string& donCfg,
    const std::string& eleCfg,
    const std::string& saveFolder,
    const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("[batchOfIVPoints]: Save folder not found");
    }

    Configuration config(
        cfg, 
        accCfg, 
        donCfg, 
        eleCfg
    );

    std::string file = saveFolder + "/" + fileName + ".npz";
    /* Output is only current right now (1D since only one output electrode) */
    std::vector<double> currentData(batchSize, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(batchSize)};
    /* Input 8D (8 electrode voltages with input_electrode set to ground) */
    std::vector<double> inputData(batchSize*config.nElectrodes, 0.0);
    std::vector<size_t> inputDataShape = {static_cast<size_t>(batchSize), config.nElectrodes};

    std::vector<double> mins(config.nElectrodes, minVoltage);
    std::vector<double> maxs(config.nElectrodes, maxVoltage);
    std::vector<std::vector<double>> samples = scaledLHC(
        batchSize,
        config.nElectrodes,
        mins,
        maxs,
        LHCSeed 
    );

    State equilState(config);
    KMCSimulator kmc(equilState);
    kmc.simulate(equilState, eqSteps, false, false);
    equilState.resetEventCounter();
    equilState.stateTime = 0.0;

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int ivPoint = 0; ivPoint < batchSize; ++ivPoint) {

            int threadSeed = threadID * 100000 + threadBaseSeed + ivPoint;
            setRandomSeed(threadSeed);
            std::vector<double> voltages = samples[ivPoint];
            voltages[outputIdx] = 0.0;
            
            double averagedCurrent = singleIVPoint(
                equilState,
                outputIdx,
                numOfTasks,
                simSteps,
                voltages
            );

            currentData[ivPoint] = averagedCurrent;
            for (int k = 0; k < config.nElectrodes; ++k) {
                inputData[k + ivPoint*config.nElectrodes] = voltages[k];
            }
        }
    }

    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "w");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "a");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");

    /* Save coords of equilState */
    std::vector<size_t> accCoordsShape = {static_cast<size_t>(equilState.nAcceptors), 2};
    std::vector<size_t> donCoordsShape = {static_cast<size_t>(equilState.nDonors), 2};
    cnpy::npz_save(file, "acc_xy", equilState.acceptorCoordinates.data(), accCoordsShape, "a");
    cnpy::npz_save(file, "don_xy", equilState.donorCoordinates.data(), donCoordsShape, "a");

    /* Neighbouring */
    std::vector<size_t> jaggedArrayLengthsShape = {static_cast<size_t>(equilState.jaggedArrayLengths.size())};
    std::vector<size_t> neighbourIndicesShape = {static_cast<size_t>(equilState.neighbourIndices.size())};
    cnpy::npz_save(file, "jagged_lengths", equilState.jaggedArrayLengths.data(), jaggedArrayLengthsShape, "a");
    cnpy::npz_save(file, "neighbour_indices", equilState.neighbourIndices.data(), neighbourIndicesShape, "a");

    /* Initial_site_energies */
    std::vector<size_t> initialSiteEnergiesShape = {static_cast<size_t>(equilState.numOfSites)};
    cnpy::npz_save(file, "init_energies", equilState.initialSiteEnergies.data(), initialSiteEnergiesShape, "a");

    /* Constant rate part */
    std::vector<double> ratePrefactors(equilState.numOfSites*equilState.numOfSites);
    std::vector<size_t> ratePrefactorsShape = {static_cast<size_t>(equilState.numOfSites), static_cast<size_t>(equilState.numOfSites)};
}

int argParser(int argc, char* argv[]) {

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
            ("accCfg", boost::program_options::value<std::string>()->required())
            ("donCfg", boost::program_options::value<std::string>()->required())
            ("eleCfg", boost::program_options::value<std::string>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())
            ("eqSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simSteps", boost::program_options::value<int>()->required())
            ("seed", boost::program_options::value<int>()->required())
            ("fileName", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
        ;

        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        std::vector<double> voltages(8, 0.0);
        if (vm.count("c_v")) {
            for (auto &s : vm["c_v"].as<std::vector<std::string>>()) {
                auto eq = s.find('=');
                int idx = std::stoi(s.substr(0, eq));
                double v = std::stod(s.substr(eq+1));
                voltages[idx] = v;
            }
        }

        singleRun(
            vm["fileName"].as<std::string>(),
            vm["eqSteps"].as<int>(),
            vm["simSteps"].as<int>(),
            voltages,
            vm["cfg"].as<std::string>(),
            vm["accCfg"].as<std::string>(),
            vm["donCfg"].as<std::string>(),
            vm["eleCfg"].as<std::string>(),
            vm["seed"].as<int>(),
            vm["saveFolder"].as<std::string>()
        );

        return 1;
    }

    if (firstCommand == "singleCurve") {

        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("numOfPoints", boost::program_options::value<int>()->required())
            ("inputIdx", boost::program_options::value<int>()->required())
            ("outputIdx", boost::program_options::value<int>()->required())
            ("minVoltage", boost::program_options::value<double>()->default_value(-1.5))
            ("maxVoltage", boost::program_options::value<double>()->default_value(1.5))
            ("eqSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simSteps", boost::program_options::value<int>()->required())
            ("numIntervals", boost::program_options::value<int>()->default_value(100))
            ("seed", boost::program_options::value<int>()->default_value(64))
            ("cfg", boost::program_options::value<std::string>()->required())
            ("accCfg", boost::program_options::value<std::string>()->required())
            ("donCfg", boost::program_options::value<std::string>()->required())
            ("eleCfg", boost::program_options::value<std::string>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())
            ("fileName", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
        ;

        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);
        
        std::vector<double> voltages(8, 0.0);
        if (vm.count("c_v")) {
            for (auto &s : vm["c_v"].as<std::vector<std::string>>()) {
                auto eq = s.find('=');
                int idx = std::stoi(s.substr(0, eq));
                double v = std::stod(s.substr(eq+1));
                voltages[idx] = v;
            }
        }

        voltages[vm["inputIdx"].as<int>()] = 0.0;
        voltages[vm["outputIdx"].as<int>()] = 0.0;  
        //std::cout << vm["cfg"].as<std::string>() << "\n";
        singleIVCurve(
            vm["numOfPoints"].as<int>(),
            vm["inputIdx"].as<int>(),
            vm["outputIdx"].as<int>(),
            vm["minVoltage"].as<double>(),
            vm["maxVoltage"].as<double>(),
            vm["eqSteps"].as<int>(),
            vm["simSteps"].as<int>(),
            vm["numIntervals"].as<int>(),
            vm["seed"].as<int>(),
            voltages,
            vm["cfg"].as<std::string>(),
            vm["accCfg"].as<std::string>(),
            vm["donCfg"].as<std::string>(),
            vm["eleCfg"].as<std::string>(),
            vm["saveFolder"].as<std::string>(),
            vm["fileName"].as<std::string>()
        );

        return 1;
    }

    if (firstCommand == "createBatch") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("batchSize", boost::program_options::value<int>()->required())
            ("minVoltage", boost::program_options::value<double>()->required())
            ("maxVoltage", boost::program_options::value<double>()->required())
            ("inputIdx", boost::program_options::value<int>()->required())
            ("outputIdx", boost::program_options::value<int>()->required())
            ("eqSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simSteps", boost::program_options::value<int>()->required())
            ("numOfTasks", boost::program_options::value<int>()->default_value(100))
            ("LHCSeed", boost::program_options::value<int>()->required())
            ("threadBaseSeed", boost::program_options::value<int>()->required())
            ("cfg", boost::program_options::value<std::string>()->required())
            ("accCfg", boost::program_options::value<std::string>()->required())
            ("donCfg", boost::program_options::value<std::string>()->required())
            ("eleCfg", boost::program_options::value<std::string>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())        
            ("fileName", boost::program_options::value<std::string>()->required())
        ;
        
        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        batchOfIVPoints(
            vm["batchSize"].as<int>(),
            vm["minVoltage"].as<double>(),
            vm["maxVoltage"].as<double>(),
            vm["inputIdx"].as<int>(),
            vm["outputIdx"].as<int>(),
            vm["eqSteps"].as<int>(),
            vm["simSteps"].as<int>(),
            vm["numOfTasks"].as<int>(),
            vm["LHCSeed"].as<int>(),
            vm["threadBaseSeed"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["accCfg"].as<std::string>(),
            vm["donCfg"].as<std::string>(),
            vm["eleCfg"].as<std::string>(),
            vm["saveFolder"].as<std::string>(),
            vm["fileName"].as<std::string>()
        );

        return 1;
    }

    return 1;
}