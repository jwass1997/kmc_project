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
    const std::string& deviceConfigs, 
    const std::string& saveFolderPath
) {

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    Configuration config(deviceConfigs, true);
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
        netEvents += (inEvents-outEvents);

        state.resetEventCounter();

        intervalCount++;
    }

    averagedCurrent = static_cast<double>(netEvents) / totalTime;

    return averagedCurrent;
}

void createDatapoint(
    const std::string& name,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& saveFolderPath
) {

}

void singleStateBatch(
    int batchSize,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& configs,
    const std::string& save,
    const std::string& batchName
) {
    if (save.empty()) {
        throw std::invalid_argument("singleStatebatch: No such folder");
    }

    std::string fileName = save + "/batch_" + batchName + ".npz";

    const int seed0 = 1234567890;

    Configuration cfg(configs, true);
    int femResolution = 1e5;

    std::vector<double> inputs(batchSize*cfg.nElectrodes, 0.0);
    std::vector<size_t> inputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(cfg.nElectrodes)};

    std::vector<double> outputs(batchSize, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize)};

    std::vector<double> voltages(batchSize, 0.0);   

    double range = maxVoltage - minVoltage;
    double vStep = range / (static_cast<double>(batchSize - 1));void singleStateBatch(
    int batchSize,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& configs,
    const std::string& save,
    const std::string& batchName
) {
    if (save.empty()) {
        throw std::invalid_argument("singleStatebatch: No such folder");
    }

    std::string fileName = save + "/batch_" + batchName + ".npz";

    const int seed0 = 1234567890;

    Configuration cfg(configs, true);
    int femResolution = 1e5;

    std::vector<double> inputs(batchSize*cfg.nElectrodes, 0.0);
    std::vector<size_t> inputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(cfg.nElectrodes)};

    std::vector<double> outputs(batchSize, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize)};

    std::vector<double> voltages(batchSize, 0.0);   

    double range = maxVoltage - minVoltage;
    double vStep = range / (batchSize - 1);

    for (int i = 0; i < batchSize; ++i) {
        voltages[i] = minVoltage + i*vStep;
    }

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed0 + static_cast<long int>(now_ns) + threadID);
        
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {

            FiniteElementeCircle fem(cfg.radius, femResolution);
            State state(cfg, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(cfg.nElectrodes, 0.void singleStateBatch(
    int batchSize,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& configs,
    const std::string& save,
    const std::string& batchName
) {
    if (save.empty()) {
        throw std::invalid_argument("singleStatebatch: No such folder");
    }

    std::string fileName = save + "/batch_" + batchName + ".npz";

    const int seed0 = 1234567890;

    Configuration cfg(configs, true);
    int femResolution = 1e5;

    std::vector<double> inputs(batchSize*cfg.nElectrodes, 0.0);
    std::vector<size_t> inputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(cfg.nElectrodes)};

    std::vector<double> outputs(batchSize, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize)};

    std::vector<double> voltages(batchSize, 0.0);   

    double range = maxVoltage - minVoltage;
    double vStep = range / (batchSize - 1);

    for (int i = 0; i < batchSize; ++i) {
        voltages[i] = minVoltage + i*vStep;
    }

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed0 + static_cast<long int>(now_ns) + threadID);
        
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {void singleStateBatch(
    int batchSize,
    int electrodeIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& configs,
    const std::string& save,
    const std::string& batchName
) {
    if (save.empty()) {
        throw std::invalid_argument("singleStatebatch: No such folder");
    }

    std::string fileName = save + "/batch_" + batchName + ".npz";

    const int seed0 = 1234567890;

    Configuration cfg(configs, true);
    int femResolution = 1e5;

    std::vector<double> inputs(batchSize*cfg.nElectrodes, 0.0);
    std::vector<size_t> inputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(cfg.nElectrodes)};

    std::vector<double> outputs(batchSize, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize)};

    double vStep = (maxVoltage - minVoltage) / (batchSize-1);

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed0 + static_cast<long int>(now_ns) + threadID);
        
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {

            FiniteElementeCircle fem(cfg.radius, femResolution);
            State state(cfg, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(cfg.nElectrodes, 0.0);
            newBoundaries[(electrodeIdx+1) % cfg.nElectrodes] = minVoltage + _batch*vStep;//minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            /* for (int i = 0; i < state.nElectrodes; ++i) {
                newBoundaries[i] = minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            } */

            state.updateBoundaries(newBoundaries, fem);

            double averagedCurrent = calculateCurrent(
                state,
                kmc,
                electrodeIdx,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals
            );

            for (int i = 0; i < cfg.nElectrodes; ++i) {
                inputs[_batch*cfg.nElectrodes + i] = newBoundaries[i];
            }

            outputs[_batch] = averagedCurrent;

            std::cout << "Finished batch#" << _batch << "\n";
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchName, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), inputShape, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), outputShape, "a");
}

            FiniteElementeCircle fem(cfg.radius, femResolution);
            State state(cfg, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(cfg.nElectrodes, 0.0);
            newBoundaries[(electrodeIdx+1) % cfg.nElectrodes] = voltages[_batch];//minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            /* for (int i = 0; i < state.nElectrodes; ++i) {
                newBoundaries[i] = minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            } */

            state.updateBoundaries(newBoundaries, fem);

            double averagedCurrent = calculateCurrent(
                state,
                kmc,
                electrodeIdx,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals
            );

            for (int i = 0; i < cfg.nElectrodes; ++i) {
                inputs[_batch*cfg.nElectrodes + i] = newBoundaries[i];
            }

            outputs[_batch] = averagedCurrent;

            std::cout << "Finished batch#" << _batch << "\n";
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchName, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), inputShape, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), outputShape, "a");
}0);
            newBoundaries[(electrodeIdx+1) % cfg.nElectrodes] = voltages[_batch];//minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            /* for (int i = 0; i < state.nElectrodes; ++i) {
                newBoundaries[i] = minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            } */

            state.updateBoundaries(newBoundaries, fem);

            double averagedCurrent = calculateCurrent(
                state,
                kmc,
                electrodeIdx,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals
            );

            for (int i = 0; i < cfg.nElectrodes; ++i) {
                inputs[_batch*cfg.nElectrodes + i] = newBoundaries[i];
            }

            outputs[_batch] = averagedCurrent;

            std::cout << "Finished batch#" << _batch << "\n";
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchName, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), inputShape, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), outputShape, "a");
}

    for (int i = 0; i < batchSize; ++i) {
        voltages[i] = minVoltage + i*vStep;
    }

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed0 + static_cast<long int>(now_ns) + threadID);
        
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {

            FiniteElementeCircle fem(cfg.radius, femResolution);
            State state(cfg, fem);
            KMCSimulator kmc(state);

            std::vector<double> newBoundaries(cfg.nElectrodes, 0.0);
            newBoundaries[(electrodeIdx+1) % cfg.nElectrodes] = voltages[_batch];//minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            /* for (int i = 0; i < state.nElectrodes; ++i) {
                newBoundaries[i] = minVoltage + (maxVoltage - minVoltage)*randomDouble01();
            } */

            state.updateBoundaries(newBoundaries, fem);

            double averagedCurrent = calculateCurrent(
                state,
                kmc,
                electrodeIdx,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals
            );

            for (int i = 0; i < cfg.nElectrodes; ++i) {
                inputs[_batch*cfg.nElectrodes + i] = newBoundaries[i];
            }

            outputs[_batch] = averagedCurrent;

            std::cout << "Finished batch#" << _batch << "\n";
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchName, {1}, "w");
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
            ("configs", boost::program_options::value<std::string>()->default_value("../configs"))
            ("save_path", boost::program_options::value<std::string>()->default_value("../data"))
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
            vm["configs"].as<std::string>(),
            vm["save_path"].as<std::string>()
        );

        return 1;
    }

    if (firstCommand == "batchRun") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("configs", boost::program_options::value<std::string>()->default_value("../configs"))
            ("save_path", boost::program_options::value<std::string>()->default_value("../data"))
            ("batchSize", boost::program_options::value<int>()->required())
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("batchName", boost::program_options::value<std::string>()->required())
        ;

        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        singleStateBatch(
            vm["batchSize"].as<int>(),
            0,
            -1.5,
            1.5,
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            100,
            vm["configs"].as<std::string>(),
            vm["save_path"].as<std::string>(),
            vm["batchName"].as<std::string>()
        );

        return 1;
    }

    return 1;
}