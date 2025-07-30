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
        netEvents += inEvents-outEvents;

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

    newBoundaries[outputIdx] = 1.2;
    for (int _v = 0; _v < numOfPoints; ++_v) {

        newBoundaries[inputIdx] = vMin + _v*vStep;
        state.updateBoundaries(newBoundaries, fem);
        
        double averageCurrent = calculateCurrent(
            state,
            kmc,
            outputIdx,
            equilibriumSteps,
            simulationSteps,
            100
        );

        outputCurrents[_v] = averageCurrent;
    }

    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "outputCurrent", outputCurrents.data(), outputShape, "a");
}

void calculateIVCurve(
    const std::string& fileName,
    int sampleSize,
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
        throw std::invalid_argument("calculateIVCurve(): Save folder not found");
    }    

    std::string filePath = saveFolderPath + "/" + fileName + ".npz";

    std::vector<double> outputs(numOfPoints*sampleSize, 0.0);

    int femRes = 1e5;
    int numOfElectrodes = 8;

    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    #pragma omp parallel
    {
        //int threadID = omp_get_thread_num();
        //auto now = std::chrono::high_resolution_clock::now();
        //auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

        #pragma omp for
        for (int s = 0; s < sampleSize; ++s) {
            setRandomSeed(seed + s); 

            Configuration config(
                cfg,
                acceptorCfg,
                donorCfg,
                electrodeCfg,
                false
            );
            FiniteElementeCircle fem(config.radius, femRes);
            State state(config, fem);
            KMCSimulator kmc(state);
            std::cout << config.energyDisorder << "\n";
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

    std::vector<double> averagedCurve(numOfPoints, 0.0);
    for (int i = 0; i < numOfPoints; ++i) {

        for (int j = 0; j < sampleSize; ++j) {
            averagedCurve[i] += outputs[i + j*numOfPoints] / sampleSize;            
        }
    }

    std::vector<size_t> shape = {static_cast<size_t>(numOfPoints)};
    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "outputs", averagedCurve.data(), shape, "a");
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
                    outputIdx,
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
    int numOfCurves,
    int inputIdx,
    int outputIdx,
    double minVoltage,
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    int seed,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolder,
    int batchID
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("singleStatebatch: Save folder not found");
    }
    std::string fileName = saveFolder + "/batch_" + std::to_string(batchID) + ".npz";

    int femResolution = 1e5;
    int numOfElectrodes = 8;

    std::vector<double> outputs(batchSize*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(numOfPoints)}; 

    double range = maxVoltage - minVoltage;
    double vStep = range / (numOfPoints - 1);

    /* std::vector<int> controlElectrodeIdx;
    for (int i = 0; i < numOfElectrodes; ++i) {
        if (i == inputIdx || i == outputIdx) {
            continue;
        }
        else {
            controlElectrodeIdx.push_back(i);
        }
    } */

    std::vector<double> mins = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
    std::vector<double> maxs = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    std::vector<std::vector<double>> voltageSample = scaledLHC(
        batchSize,
        numOfElectrodes,
        mins,
        maxs            
    );

    #pragma omp parallel
    {   
        #pragma omp for
        for (int _batch = 0; _batch < batchSize; ++_batch) {

            int threadID = omp_get_thread_num();
            //auto now = std::chrono::high_resolution_clock::now();
            //auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

            for (int _c = 0; _c < numOfCurves; ++_c) {

                setRandomSeed(seed + _batch * numOfCurves + _c);

                Configuration config(
                    cfg, 
                    acceptorCfg, 
                    donorCfg, 
                    electrodeCfg, 
                    false
                );                

                FiniteElementeCircle fem(config.radius, femResolution);
                State state(config, fem);
                KMCSimulator kmc(state);

                std::vector<double> newBoundaries(numOfElectrodes, 0.0);
                for (int i = 0; i < numOfElectrodes; ++i) {
                    newBoundaries[i] = voltageSample[_batch][i];
                }
                newBoundaries[outputIdx] = 0.0; 
                for (int _v = 0; _v < numOfPoints; _v++) {

                    newBoundaries[inputIdx] = minVoltage + _v*vStep;
                    state.updateBoundaries(newBoundaries, fem);

                    double averagedCurrent = calculateCurrent(
                        state,
                        kmc,
                        outputIdx,
                        equilibriumSteps,
                        simulationSteps,
                        numOfIntervals
                    );

                    outputs[_batch*numOfPoints + _v] += averagedCurrent / static_cast<double>(numOfCurves);
                }
            }
        }
    }
    
    cnpy::npz_save(fileName, "ID", &batchID, {1}, "w");
    cnpy::npz_save(fileName, "out", outputs.data(), outputShape, "a");
}

void lineSweep(
    const std::string& fileName,
    const std::string& constParamName,
    const std::string& varParamName,
    double constParam,
    double varParamMin,
    double varParamMax,
    int N,
    int sampleSize,
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
        throw std::invalid_argument("lineSweep(): No save folder found");
    }
    std::string filePath = saveFolderPath + "/" + fileName + ".npz";

    int femRes = 1e5;
    int numOfElectrodes = 8;

    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    double paramStep = (varParamMax - varParamMin) / static_cast<double>(N-1);

    std::vector<double> outputs(N*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(N), static_cast<size_t>(numOfPoints)};

    #pragma omp parallel
    {   
        #pragma omp for
        for (int i = 0; i < N; ++i) {
            
            for (int s = 0; s < sampleSize; ++s) {

                setRandomSeed(seed + i * sampleSize + s);

                Configuration config(
                    cfg,
                    acceptorCfg,
                    donorCfg,
                    electrodeCfg,
                    false
                );

                if (constParamName == "sigma") {
                    config.energyDisorder = constParam*e / config.kbT;                
                }
                else if(constParamName == "temp") {
                    config.T = constParam;
                }

                double varParamVal = varParamMin + i*paramStep;
                if (varParamName == "sigma") {
                    config.energyDisorder = varParamVal*e / config.kbT;
                }
                else if ( varParamName == "temp") {
                    config.T = varParamVal;
                }

                FiniteElementeCircle fem(config.radius, femRes);
                State state(config, fem);
                KMCSimulator kmc(state);

                std::vector<double> newBoundaries(numOfElectrodes, 0.0);
                for (int _b = 0; _b < newBoundaries.size(); ++_b) {
                    newBoundaries[_b] = voltages[_b];
                }
                newBoundaries[outputIdx] = 0.0;

                for (int _v = 0; _v < numOfPoints; _v++) {

                    newBoundaries[inputIdx] = vMin + _v*vStep;
                    state.updateBoundaries(newBoundaries, fem);

                    double current = calculateCurrent(
                        state,
                        kmc,
                        outputIdx,
                        equilibriumSteps,
                        simulationSteps,
                        100
                    );

                    outputs[i*numOfPoints + _v] += current / static_cast<double>(sampleSize);
                }
            }
        }
    }
    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "out", outputs.data(), outputShape, "a");
}

void param2DSweep(
    const std::string& fileName,
    const std::string& paramName1,
    const std::string& paramName2,
    double paramValue1,
    double paramValue2,
    int sampleSize,
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
        throw std::invalid_argument("param2DLineSweep(): No save folder found");
    }
    std::string filePath = saveFolderPath + "/" + fileName + ".npz";

    std::vector<double> outputs(numOfPoints*sampleSize, 0.0);

    int femRes = 1e5;
    int numOfElectrodes = 8;

    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        //auto now = std::chrono::high_resolution_clock::now();
        //auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(seed + threadID); 

        #pragma omp for
        for (int s = 0; s < sampleSize; ++s) {

            Configuration config(
                cfg,
                acceptorCfg,
                donorCfg,
                electrodeCfg,
                false
            );
            
            if (paramName1 == "sigma") {
                config.energyDisorder = paramValue1*e / config.kbT;                
            }
            else if(paramName1 == "temp") {
                config.T = paramValue1;
            }

            if (paramName2 == "temp") {
                config.T = paramValue2;                
            }
            else if(paramName2 == "sigma") {
                config.energyDisorder = paramValue2*e / config.kbT;
            }

            FiniteElementeCircle fem(config.radius, femRes);
            State state(config, fem);
            KMCSimulator kmc(state);
            //std::cout << config.energyDisorder << "\n";
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
                    outputIdx,
                    equilibriumSteps,
                    simulationSteps,
                    100
                );

                outputs[s*numOfPoints + _v] = averagedCurrent;
            }
        }
    }

    std::vector<double> averagedCurve(numOfPoints, 0.0);
    for (int i = 0; i < numOfPoints; ++i) {

        for (int j = 0; j < sampleSize; ++j) {
            averagedCurve[i] += outputs[i + j*numOfPoints] / sampleSize;            
        }
    }

    std::vector<size_t> shape = {static_cast<size_t>(numOfPoints)};
    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "outputs", averagedCurve.data(), shape, "a");
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
            ("savePath", boost::program_options::value<std::string>()->required())
            ("batchSize", boost::program_options::value<int>()->required())
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("numOfCurves", boost::program_options::value<int>()->required())
            ("inputIdx", boost::program_options::value<int>()->required())
            ("outputIdx", boost::program_options::value<int>()->required())
            ("seed", boost::program_options::value<int>()->required())
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
            vm["numOfCurves"].as<int>(),
            vm["inputIdx"].as<int>(),
            vm["outputIdx"].as<int>(),
            -1.5,
            1.5,
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            100,
            vm["seed"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["savePath"].as<std::string>(),
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

    if (firstCommand == "calcIVCurve") {
        
        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("file_name", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
            ("sampleSize", boost::program_options::value<int>()->required())
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
        
        calculateIVCurve(
            vm["file_name"].as<std::string>(),
            vm["sampleSize"].as<int>(),
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

    if (firstCommand == "2DSweep") {
        
        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("file_name", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
            ("paramName1", boost::program_options::value<std::string>()->required())
            ("paramName2", boost::program_options::value<std::string>()->required())
            ("paramValue1", boost::program_options::value<double>()->required())
            ("paramValue2", boost::program_options::value<double>()->required())
            ("sampleSize", boost::program_options::value<int>()->required())
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
        
        param2DSweep(
            vm["file_name"].as<std::string>(),
            vm["paramName1"].as<std::string>(),
            vm["paramName2"].as<std::string>(),
            vm["paramValue1"].as<double>(),
            vm["paramValue2"].as<double>(),
            vm["sampleSize"].as<int>(),
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

    if (firstCommand == "lineSweep") {
        
        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("file_name", boost::program_options::value<std::string>()->required())
            ("c_v", boost::program_options::value<std::vector<std::string>>()->composing(), "electrode index=value")
            ("constParamName", boost::program_options::value<std::string>()->required())
            ("varParamName", boost::program_options::value<std::string>()->required())
            ("constParam", boost::program_options::value<double>()->required())
            ("varParamMin", boost::program_options::value<double>()->required())
            ("varParamMax", boost::program_options::value<double>()->required())
            ("N", boost::program_options::value<int>()->required())
            ("sampleSize", boost::program_options::value<int>()->required())
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
        
        lineSweep(
            vm["file_name"].as<std::string>(),
            vm["constParamName"].as<std::string>(),
            vm["varParamName"].as<std::string>(),
            vm["constParam"].as<double>(),
            vm["varParamMin"].as<double>(),
            vm["varParamMax"].as<double>(),
            vm["N"].as<int>(),
            vm["sampleSize"].as<int>(),
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