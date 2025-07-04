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

std::vector<double> calculateSlopes(
    std::vector<double> fX,
    std::vector<double> X
) {
        
    if (fX.size() != X.size()) {
        std::cerr << "Slope can not be calculated if sizes do not match" << "\n";
    }

    std::vector<double> slopes(fX.size(), 0.0);

    double _slope = 0.0;
    for (int i = 1; i < slopes.size(); ++i) {
        _slope = (fX[i] - fX[i-1]) / (X[i] - X[i-1]);
        slopes[i] = _slope;
    }

    return slopes;
}

void createDirectoryFromStringPath(const std::string& path, const std::string& folderName) {

    if (folderName.empty()) {
        std::cerr << "Must specify a folder name " << "\n";
        return;
    }

    std::filesystem::path directoryPath = path.empty() ? std::filesystem::current_path() : std::filesystem::path(path);
    std::filesystem::path newFolder = std::filesystem::path(directoryPath)/folderName;

    try {
        if(std::filesystem::create_directory(newFolder)) {
            std::cout << "Folder has been created: " << newFolder << "\n";
        }
        else {
            std::cout << "Folder already exists: " << newFolder << "\n";
        }
    }
    catch(const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating new folder: " << e.what() << "\n";
    }
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
    FiniteElementeCircle fem(config.radius, 1e5);

    State state(config, fem);
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

double calculateCurrent(
    State& state,
    KMCSimulator& kmc,
    int electrodeIdx,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals
) {

    double averagedCurrent = 0.0;
    double totalTime = 0.0;
    int netEvents = 0;
    int intervalSteps = simulationSteps / numOfIntervals;
    int intervalCount = 0;

    state.resetEventCounter();
    kmc.simulate(state, equilibriumSteps, false, false);

    while (intervalCount < numOfIntervals) {
        //std::cout << "Thread #" << omp_get_thread_num() << " at count" << intervalCount << "\n";
        double startClock = state.stateTime;
        kmc.simulate(state, intervalSteps, false, true);
        double endClock = state.stateTime;

        double elapsedTime = endClock - startClock;
        int inEvents = 0;
        int outEvents = 0;
        for (int i = 0; i < state.numOfSites; ++i) {
            outEvents += state.eventCounter[(electrodeIdx+state.nAcceptors)*state.numOfSites + i];
            inEvents += state.eventCounter[state.numOfSites*i + (electrodeIdx+state.nAcceptors)];
        }
        totalTime += elapsedTime;
        netEvents += -(inEvents-outEvents);

        state.resetEventCounter();

        intervalCount++;
    }

    averagedCurrent = static_cast<double>(netEvents) / totalTime;

    return averagedCurrent;
}

void oneDimensionalCurve(
    const std::string& fileName,
    std::vector<double> voltages,
    int inputIdx,
    int outputIdx,
    double vMin,
    double vMax,
    int numOfPoints,
    int equilibriumSteps,
    int simulationSteps,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolderPath,
    bool randomGeometry,
    int seed
) {

    if (saveFolderPath.empty()) {
        throw std::invalid_argument("createDatapoint(): No such folder");
    }

    std::string filePath = saveFolderPath + "/data_point_" + fileName + ".npz";
    
    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    std::vector<double> outputCurrents(numOfPoints, 0.0);
    std::vector<size_t> outputShape{static_cast<size_t>(numOfPoints)};

    std::vector<double> newBoundaries(8, 0.0);
    for (int i = 0; i < newBoundaries.size(); ++i) {
        newBoundaries[i] = voltages[i];
        std::cout << voltages[i] << "\n";
    }

    setRandomSeed(seed);

    int femRes = 1e5;
    Configuration config(cfg, acceptorCfg, donorCfg, electrodeCfg, randomGeometry);
    FiniteElementeCircle fem(config.radius, femRes);
    State state(config, fem);
    KMCSimulator kmc(state);

    newBoundaries[outputIdx] = 0.0;
    for (int _v = 0; _v < numOfPoints; ++_v) {

        newBoundaries[inputIdx] = vMin + _v*vStep;
        state.updateBoundaries(newBoundaries, fem);
        
        double averageCurrent = calculateCurrent(
            state,
            kmc,
            inputIdx,
            equilibriumSteps,
            simulationSteps,
            100
        );

        outputCurrents[_v] = averageCurrent;
    }

    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "outputCurrent", outputCurrents.data(), outputShape, "a");
}

void sampleOfIVCurves(
    int trial,
    std::vector<double> voltages,
    int numOfSamples,
    int inputIdx,
    int outputIdx,
    double vMin,
    double vMax,
    int numOfPoints,
    int equilibriumSteps,
    int simulationSteps,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolder,
    int seed
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("sampleOfIVCurves(): Save folder not found");
    }

    std::string fileName = saveFolder + "/trial_" + std::to_string(trial) + "/curves_trial=" + std::to_string(trial) + ".npz";

    const int seed0 = 1234567890;

    int femRes = 1e5;
    int numOfElectrodes = 8;

    std::vector<double> outputs(numOfSamples*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(numOfSamples), static_cast<size_t>(numOfPoints)};

    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed + threadID);

        #pragma omp for
        for (int s = 0; s < numOfSamples; ++s) {

            auto acceptorDist = acceptorCfg + "/acceptors_trial=" + std::to_string(trial) + "_sample=" + std::to_string(s);

            Configuration config(cfg, acceptorDist, donorCfg, electrodeCfg, false);
            FiniteElementeCircle fem(config.radius, femRes);
            State state(config, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(numOfElectrodes, 0.0);
            for (int i = 0; i < newBoundaries.size(); ++i) {
                newBoundaries[i] = voltages[i];
            }

            newBoundaries[outputIdx] = 0.0;
            for (int _v = 0; _v < numOfPoints; _v++) {

                newBoundaries[inputIdx] = vMin + _v*vStep;
                state.updateBoundaries(newBoundaries, fem);

                double averagedCurrent = calculateCurrent(
                    state,
                    kmc,
                    inputIdx,
                    equilibriumSteps,
                    simulationSteps,
                    100
                );

                outputs[s*numOfPoints + _v] = averagedCurrent;
            }
        }
    }
    cnpy::npz_save(fileName, "ID", &fileName, {1}, "w");
    cnpy::npz_save(fileName, "outputs", outputs.data(), outputShape, "a");
}

void batchOfIVCurves(
    int batchSize,
    int numOfPoints,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& save,
    int batchID
) {
    if (save.empty()) {
        throw std::invalid_argument("singleStatebatch: Save folder not found");
    }

    std::string fileName = save + "/batch_" + std::to_string(batchID) + ".npz";

    const int seed0 = 1234567890;

    int femResolution = 1e5;
    int numOfElectrodes = 8;

    std::vector<double> inputs(batchSize*numOfElectrodes, 0.0);
    std::vector<size_t> inputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(numOfElectrodes)};

    std::vector<double> outputs(batchSize*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(numOfPoints)}; 

    double range = maxVoltage - minVoltage;
    double vStep = range / (numOfPoints - 1);

    #pragma omp parallel
    {   

        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed0 + static_cast<long int>(now_ns) + threadID);

        Configuration config(cfg, acceptorCfg, donorCfg, electrodeCfg, false);
        
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {

            
            FiniteElementeCircle fem(config.radius, femResolution);
            State state(config, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(numOfElectrodes, 0.0);
            for (int i = 0; i < numOfElectrodes; ++i) {
                inputs[_batch*numOfElectrodes + i] = minVoltage + (maxVoltage - minVoltage)*randomDouble01();//newBoundaries[i];    
            }

            for (int _v = 0; _v < numOfPoints; _v++) {

                newBoundaries[(electrodeIdx+1) % numOfElectrodes] = minVoltage + _v*vStep;//minVoltage + (maxVoltage - minVoltage)*randomDouble01();
                state.updateBoundaries(newBoundaries, fem);

                double averagedCurrent = calculateCurrent(
                    state,
                    kmc,
                    electrodeIdx,
                    equilibriumSteps,
                    simulationSteps,
                    numOfIntervals
                );

                outputs[_batch*numOfPoints + _v] = averagedCurrent;
                //std::cout << minVoltage + _batch*vStep << "\n";
            }
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchID, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), inputShape, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), outputShape, "a");
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
                  << " singleRun --configs <string> --save_path <string> --equilibriumSteps <int> --simulationSteps <int> --deviceName <string>\n"
                  << " batchRun --configs <string> --save_path <string> --batchSize <int> --equilibriumSteps <int> --simulationSteps <int> --batchName <string>\n";
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

    if (firstCommand == "batchRun") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("cfg", boost::program_options::value<std::string>()->required())
            ("acceptorCfg", boost::program_options::value<std::string>()->required())
            ("donorCfg", boost::program_options::value<std::string>()->required())
            ("electrodeCfg", boost::program_options::value<std::string>()->required())
            ("save_path", boost::program_options::value<std::string>()->required())
            ("batchSize", boost::program_options::value<int>()->required())
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("batchID", boost::program_options::value<int>()->required())
        ;
        
        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        batchOfIVCurves(
            vm["batchSize"].as<int>(),
            vm["numOfPoints"].as<int>(),
            0,
            -1.5,
            1.5,
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            100,
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["save_path"].as<std::string>(),
            vm["batchID"].as<int>()
        );

        return 1;
    }
    
    if (firstCommand == "findControlVoltages") {
        
        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("file_name", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
            ("input_idx", boost::program_options::value<int>()->required())
            ("output_idx", boost::program_options::value<int>()->required())
            ("vMin", boost::program_options::value<double>()->default_value(-1.5))
            ("vMax", boost::program_options::value<double>()->default_value(1.5))
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("cfg", boost::program_options::value<std::string>()->required())
            ("acceptorCfg", boost::program_options::value<std::string>()->required())
            ("donorCfg", boost::program_options::value<std::string>()->required())
            ("electrodeCfg", boost::program_options::value<std::string>()->required())
            ("saveFolderPath", boost::program_options::value<std::string>()->required())
            ("randomGeometry", boost::program_options::value<int>()->default_value(0))
            ("seed", boost::program_options::value<int>()->default_value(64))
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

        voltages[vm["input_idx"].as<int>()] = 0.0;
        voltages[vm["output_idx"].as<int>()] = 0.0;  
        
        oneDimensionalCurve(
            vm["file_name"].as<std::string>(),
            voltages,
            vm["input_idx"].as<int>(),
            vm["output_idx"].as<int>(),
            vm["vMin"].as<double>(),
            vm["vMax"].as<double>(),
            vm["numOfPoints"].as<int>(),
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["saveFolderPath"].as<std::string>(),
            vm["randomGeometry"].as<int>(),
            vm["seed"].as<int>()
        );

        return 1;
    }

    if (firstCommand == "findNodeDistribution") {
        
        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("numOfSamples", boost::program_options::value<int>()->required())
            ("trial", boost::program_options::value<int>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
            ("input_idx", boost::program_options::value<int>()->required())
            ("output_idx", boost::program_options::value<int>()->required())
            ("vMin", boost::program_options::value<double>()->default_value(-1.5))
            ("vMax", boost::program_options::value<double>()->default_value(1.5))
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("cfg", boost::program_options::value<std::string>()->required())
            ("acceptorCfg", boost::program_options::value<std::string>()->required())
            ("donorCfg", boost::program_options::value<std::string>()->required())
            ("electrodeCfg", boost::program_options::value<std::string>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())
            ("randomGeometry", boost::program_options::value<int>()->default_value(0))
            ("seed", boost::program_options::value<int>()->default_value(64))
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

        voltages[vm["input_idx"].as<int>()] = 0.0;
        voltages[vm["output_idx"].as<int>()] = 0.0;  
        
        sampleOfIVCurves(
            vm["trial"].as<int>(),
            voltages,
            vm["numOfSamples"].as<int>(),
            vm["input_idx"].as<int>(),
            vm["output_idx"].as<int>(),
            vm["vMin"].as<double>(),
            vm["vMax"].as<double>(),
            vm["numOfPoints"].as<int>(),
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["saveFolder"].as<std::string>(),
            vm["seed"].as<int>()
        );

        return 1;
    }
    
    return 1;
}