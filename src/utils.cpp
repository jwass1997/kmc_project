#include <cstddef>
#include <random>
#include <algorithm>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>

#include "utils.h"
#include "State.h"
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

uint64_t mix(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
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

    Configuration config(
        cfg, 
        acceptorCfg, 
        donorCfg, 
        electrodeCfg
    );

    State state(config, seed);
    KMCSimulator kmc(state, seed);
    kmc.simulate(state, eqSteps, false, false);
    state.resetEventCounter();
    state.stateTime = 0.0;

    state.updateBoundaries(voltages);
    kmc.simulate(state, simSteps, false, true);

    int nAcceptors = state.nAcceptors;
    int nElectrodes = state.nElectrodes;
    int nDonors = state.nDonors;

    std::vector<double> flattenedAcceptorCoordinates(2*state.nAcceptors, 0.0);
    std::vector<double> flattenedDonorCoordinates(2*state.nDonors, 0.0);
    std::vector<double> flattenedElectrodeCoordinates(2*state.nElectrodes, 0.0);    
    std::vector<int> flattenedEventCounts(state.numOfSites*state.numOfSites, 0);
    std::vector<double> flattenedEnergies(state.numOfSites, 0.0);

    std::vector<size_t> shapeFlattenedAcceptorCoordinates = {static_cast<size_t>(nAcceptors), 2};
    std::vector<size_t> shapeFlattenedDonorCoordinates = {static_cast<size_t>(nDonors), 2};
    std::vector<size_t> shapeFlattenedElectrodeCoordinates = {static_cast<size_t>(nElectrodes), 2};
    std::vector<size_t> shapeFlattenedEventCounts = {static_cast<size_t>(nAcceptors+nElectrodes), static_cast<size_t>(nAcceptors+nElectrodes)};
    std::vector<size_t> shapeFlattenedEnergies = {static_cast<size_t>(state.numOfSites)};

    std::string deviceName = saveFolderPath + "/" + ID + ".npz";
    cnpy::npz_save(deviceName, "ID", &ID, {1}, "w"); 

    for (int i = 0; i < nAcceptors; ++i) 
    {
        flattenedAcceptorCoordinates[i*2] = state.acceptorCoordinates[i*2];
        flattenedAcceptorCoordinates[i*2 + 1] = state.acceptorCoordinates[i*2 + 1];
    }
    for (int i = 0; i < nDonors; ++i) 
    {
        flattenedDonorCoordinates[i*2] = state.donorCoordinates[i*2];
        flattenedDonorCoordinates[i*2 + 1] = state.donorCoordinates[i*2 + 1];
    }
    for (int i = 0; i < nElectrodes; ++i) 
    {
        flattenedElectrodeCoordinates[i*2] = state.electrodeCoordinates[i*2];
        flattenedElectrodeCoordinates[i*2 + 1] = state.electrodeCoordinates[i*2 + 1];
    }

    for (int j = 0; j < nAcceptors+nElectrodes; ++j) 
    {
        for (int i = 0; i <nAcceptors+nElectrodes; ++i) 
        {
            flattenedEventCounts[j*state.numOfSites + i] = state.eventCounter[j*state.numOfSites + i];
        }
    }

    for (int i = 0; i < state.numOfSites; ++i)
    {
        flattenedEnergies[i] = state.siteEnergies[i];
    }

    double total_time = state.stateTime;

    cnpy::npz_save(deviceName, "acc_xy", flattenedAcceptorCoordinates.data(), shapeFlattenedAcceptorCoordinates, "a");
    cnpy::npz_save(deviceName, "don_xy", flattenedDonorCoordinates.data(), shapeFlattenedDonorCoordinates, "a");
    cnpy::npz_save(deviceName, "ele_xy", flattenedElectrodeCoordinates.data(), shapeFlattenedElectrodeCoordinates, "a");
    cnpy::npz_save(deviceName, "event_matrix", flattenedEventCounts.data(), shapeFlattenedEventCounts, "a");
    cnpy::npz_save(deviceName, "energies", flattenedEnergies.data(), shapeFlattenedEnergies, "a");
    cnpy::npz_save(deviceName, "sim_time", &total_time, {1}, "a");
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

    uint64_t kmc_seed = mix(
        uint64_t(0x9e3779b97f4a7c15ULL) ^
        uint64_t(outputIdx) ^
        uint64_t(simSteps) ^
        uint64_t(numOfTasks)
    );

    KMCSimulator kmc(equilState, kmc_seed);

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
    int inputIdx, int outputIdx,
    double minVoltage, double maxVoltage,
    int eqSteps, int simSteps, int numOfTasks,
    int seed,
    std::vector<double> controlVoltages,
    const std::string& cfg, const std::string& accCfg, const std::string& donCfg, const std::string& eleCfg,
    const std::string& saveFolder, const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("[singleIVCurve]: Save folder not found");
    }

    std::string file = saveFolder + "/" + fileName + ".npz";

    std::vector<double> currentData(numOfPoints, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(numOfPoints)};

    std::vector<double> currentStd(numOfPoints, 0.0);
    std::vector<size_t> currentStdShape = {static_cast<size_t>(numOfPoints)};

    std::vector<double> currentSem(numOfPoints, 0.0);
    std::vector<size_t> currentSemShape = {static_cast<size_t>(numOfPoints)};

    std::vector<double> controlData(controlVoltages.size(), 0.0);
    std::vector<size_t> controlDataShape = {controlVoltages.size()};

    for (int cv = 0; cv < controlVoltages.size(); ++cv) {
        controlData[cv] = controlVoltages[cv];
    }
    controlData[inputIdx] = -999.999;
    controlData[outputIdx] = -999.999;

    double range = maxVoltage - minVoltage;
    double vStep = range / static_cast<double>(numOfPoints-1);

    Configuration config(
        cfg,
        accCfg,
        donCfg,
        eleCfg
    );

    uint64_t stateSeed = mix(
        uint64_t(0x9e3779b97f4a7c15ULL) ^
        uint64_t(seed) ^
        uint64_t(simSteps) ^
        uint64_t(numOfTasks)
    );
    State state(config, stateSeed);

    uint64_t kmc_init_seed = mix(
        uint64_t(0x9e3779b97f4a7c15ULL) ^
        uint64_t(stateSeed) ^
        uint64_t(simSteps) ^
        uint64_t(numOfTasks)
    );
    KMCSimulator kmc_init(state, kmc_init_seed);
    kmc_init.simulate(state, eqSteps, false, false);
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
            /* double currentOutput = singleIVPoint(
                state, 
                outputIdx, 
                numOfTasks, 
                simSteps, 
                voltages
            ); */

            State equilState(state);
            equilState.updateBoundaries(voltages);

            uint64_t kmc_seed = mix(
                uint64_t(0x9e3779b97f4a7c15ULL) ^
                uint64_t(outputIdx) ^
                uint64_t(simSteps) ^
                uint64_t(numOfTasks) ^
                uint64_t(threadID)
            );

            KMCSimulator kmc(equilState, kmc_seed);

            equilState.resetEventCounter();

            double averagedCurrent = 0.0;
            double totalTime = 0.0;
            int intervalSteps = simSteps / numOfTasks;
            int netEvents = 0;

            double meanW = 0.0;
            double M2w = 0.0;
            double wSum = 0.0;
            double w2Sum = 0.0;

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

                double Ii = e * static_cast<double>(inEvents-outEvents) / elapsedTime;
                averagedCurrent += Ii*elapsedTime;
                
                /**
                 * weighted Welford update
                 */
                double w = elapsedTime;
                double prevMean = meanW;
                wSum += w;
                meanW += w * (Ii - meanW) / wSum;
                M2w += w * (Ii - prevMean) * (Ii - meanW);
                w2Sum += w * w;

                equilState.resetEventCounter();

                intervalCount++;
            }

            averagedCurrent /= totalTime;

            double weightedVar = (wSum > 0.0) ? (M2w / wSum) : 0.0;
            double sampleStd = std::sqrt(weightedVar);
            double sem = (wSum > 0.0) ? (sampleStd * std::sqrt(w2Sum) / wSum) : 0.0;

            currentData[v] = averagedCurrent;
            currentStd[v] = sampleStd;
            currentSem[v] = sem;
        }
    }

    cnpy::npz_save(file, "current", currentData.data(), currentDataShape, "w");
    cnpy::npz_save(file, "currentStd", currentStd.data(), currentStdShape, "a");
    cnpy::npz_save(file, "sem", currentSem.data(), currentSemShape, "a");
    cnpy::npz_save(file, "control", controlData.data(), controlDataShape, "a");
    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "a");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
}

void batchFromSingleState(
    int batchSize,
    double minVoltage,
    double maxVoltage,
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
        throw std::invalid_argument("[batchFromSingleState]: Save folder not found");
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

    std::vector<double> currentStd(batchSize, 0.0);
    std::vector<size_t> currentStdShape = {static_cast<size_t>(batchSize)};

    std::vector<double> currentSem(batchSize, 0.0);
    std::vector<size_t> currentSemShape = {static_cast<size_t>(batchSize)};

    /* Input 8D (8 electrode voltages with input_electrode set to ground) */
    std::vector<double> inputData(batchSize*config.nElectrodes, 0.0);
    std::vector<size_t> inputDataShape = {static_cast<size_t>(batchSize), static_cast<size_t>(config.nElectrodes)};
    /* LHC sampled voltages */
    std::vector<double> mins(config.nElectrodes, minVoltage);
    std::vector<double> maxs(config.nElectrodes, maxVoltage);
    std::vector<std::vector<double>> samples = scaledLHC(
        batchSize,
        config.nElectrodes,
        mins,
        maxs,
        LHCSeed 
    );
    /* Shared initial state */
    uint64_t init_state_seed = mix(uint64_t(threadBaseSeed) ^ 0xA54FF53A5F1D36F1ULL);
    State initState(config, init_state_seed);

    uint64_t init_kmc_seed = mix(init_state_seed ^ 0x9E3779B97F4A7C15ULL);
    KMCSimulator init_kmc(initState, init_kmc_seed);

    init_kmc.simulate(initState, eqSteps, false, false);
    initState.resetEventCounter();
    initState.stateTime = 0.0;

    /* Vector for saving initial site energies */
    std::vector<double> initialEnergies(batchSize*initState.numOfSites, 0.0);

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int ivPoint = 0; ivPoint < batchSize; ++ivPoint) {

            uint64_t kmcSeed = mix(
                uint64_t(threadBaseSeed) ^
                (uint64_t(threadID) << 32) ^
                uint64_t(ivPoint) ^
                0xD1B54A32D192ED03ULL
            );

            std::vector<double> voltages = samples[ivPoint];
            voltages[outputIdx] = 0.0;
            
            State equilState(initState);

            equilState.updateBoundaries(voltages);

            for (int _s = 0; _s < equilState.numOfSites; ++_s) {
                initialEnergies[_s + ivPoint*equilState.numOfSites] = equilState.siteEnergies[_s];
            }

            KMCSimulator kmc(equilState, kmcSeed);

            equilState.resetEventCounter();

            double averagedCurrent = 0.0;
            double totalTime = 0.0;
            int intervalSteps = simSteps / numOfTasks;
            int netEvents = 0;

            double meanW = 0.0;
            double M2w = 0.0;
            double wSum = 0.0;
            double w2Sum = 0.0;

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

                double Ii = static_cast<double>(inEvents-outEvents) / elapsedTime;
                averagedCurrent += Ii*elapsedTime;
                
                /**
                 * weighted Welford update
                 */
                double w = elapsedTime;
                double prevMean = meanW;
                wSum += w;
                meanW += w * (Ii - meanW) / wSum;
                M2w += w * (Ii - prevMean) * (Ii - meanW);
                w2Sum += w * w;
                
                equilState.resetEventCounter();
                ++intervalCount;
            }

            averagedCurrent /= totalTime;
            double weightedVar = (wSum > 0.0) ? (M2w / wSum) : 0.0;
            double sampleStd = std::sqrt(weightedVar);
            double sem = (wSum > 0.0) ? (sampleStd * std::sqrt(w2Sum) / wSum) : 0.0;

            currentData[ivPoint] = averagedCurrent;
            currentStd[ivPoint] = sampleStd;
            currentSem[ivPoint] = sem;
            for (int k = 0; k < config.nElectrodes; ++k) {
                inputData[k + ivPoint*config.nElectrodes] = voltages[k];
            }
        }
    }
    /* Input-Output data */
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "w");
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "a");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");
    cnpy::npz_save(file, "sampleStd", currentStd.data(), currentStdShape, "a");
    cnpy::npz_save(file, "sem", currentSem.data(), currentSemShape, "a");

    /* Save coords of equilState */
    std::vector<size_t> accCoordsShape = {static_cast<size_t>(initState.nAcceptors), 2};
    std::vector<size_t> donCoordsShape = {static_cast<size_t>(initState.nDonors), 2};
    std::vector<size_t> eleCoordsShape = {static_cast<size_t>(initState.nElectrodes), 2};
    cnpy::npz_save(file, "acc_xy", initState.acceptorCoordinates.data(), accCoordsShape, "a");
    cnpy::npz_save(file, "don_xy", initState.donorCoordinates.data(), donCoordsShape, "a");
    cnpy::npz_save(file, "ele_xy", initState.electrodeCoordinates.data(), eleCoordsShape, "a");

    /* Neighbouring */
    std::vector<size_t> jaggedArrayLengthsShape = {static_cast<size_t>(initState.jaggedArrayLengths.size())};
    std::vector<size_t> neighbourIndicesShape = {static_cast<size_t>(initState.neighbourIndices.size())};
    cnpy::npz_save(file, "jagged_lengths", initState.jaggedArrayLengths.data(), jaggedArrayLengthsShape, "a");
    cnpy::npz_save(file, "neighbour_indices", initState.neighbourIndices.data(), neighbourIndicesShape, "a");

    /* Initial site energies and different energy parts */
    std::vector<size_t> initialEnergiesShape = {static_cast<size_t>(batchSize), static_cast<size_t>(initState.numOfSites)};
    std::vector<size_t> accDonInteractionShape = {static_cast<size_t>(initState.nAcceptors)};
    std::vector<size_t> accInteractionShape = {static_cast<size_t>(initState.nAcceptors)};
    std::vector<size_t> randomEnergiesShape = {static_cast<size_t>(initState.nAcceptors)};
    cnpy::npz_save(file, "init_energies", initialEnergies.data(), initialEnergiesShape, "a");
    cnpy::npz_save(file, "acc_don_int", initState.acceptorDonorInteraction.data(), accDonInteractionShape, "a");
    cnpy::npz_save(file, "acc_acc_int", initState.acceptorInteraction.data(), accInteractionShape, "a");
    cnpy::npz_save(file, "rand_energies", initState.randomEnergies.data(), randomEnergiesShape, "a");

    /* Constant rate part */
    std::vector<double> ratePrefactors(initState.numOfSites*initState.numOfSites, 0.0);
    std::vector<size_t> ratePrefactorsShape = {static_cast<size_t>(initState.numOfSites), static_cast<size_t>(initState.numOfSites)};
    /* Reconstructing NxN matrix from jagged arrays */
    for (int l = 0; l < initState.jaggedArrayLengths.size()-1; ++l) {

        int jaggedLength = initState.jaggedArrayLengths[l+1] - initState.jaggedArrayLengths[l];

        for (int m = initState.jaggedArrayLengths[l]; m < initState.jaggedArrayLengths[l+1]; ++m) {

            ratePrefactors[initState.neighbourIndices[m] + l*initState.numOfSites] = init_kmc.constantTransitionRates[m];

        }
    }
    cnpy::npz_save(file, "rate_prefactors", ratePrefactors.data(), ratePrefactorsShape, "a");
}

void batchOfMultipleStates(
    int batchSize,
    double minVoltage, double maxVoltage,
    int nAcceptors, int nElectrodes, int nDonors,
    double radius, double nu0, double a, double T,
    double energyDisorder,
    double electrodeWidth,
    double minHopDistance, double maxHopDistance,
    int Nr, int Nt,
    std::string distType,
    int outputIdx,
    int eqSteps, int simSteps, int numOfTasks,
    int LHCSeed, int threadBaseSeed,
    const std::string& saveFolder, const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("[batchOfMultipleStates]: Save folder not found");
    }

    int numOfSites = nAcceptors + nElectrodes;

    /* Energy contributions */
    std::vector<double> initEn(batchSize * numOfSites, 0.0);
    std::vector<size_t> initEnShape = {size_t(batchSize), size_t(numOfSites)};

    std::vector<double> randEn(batchSize * nAcceptors, 0.0);
    std::vector<size_t> randEnShape = {size_t(batchSize), size_t(nAcceptors)};

    std::vector<double> accDonInt(batchSize * nAcceptors, 0.0);
    std::vector<size_t> accDonIntShape = {size_t(batchSize), size_t(nAcceptors)};

    std::vector<double> accAccInt(batchSize * nAcceptors, 0.0);
    std::vector<size_t> accAccIntShape = {size_t(batchSize), size_t(nAcceptors)};


    std::string file = saveFolder + "/" + fileName + ".npz";
    /* Output current */
    std::vector<double> currentData(batchSize * nElectrodes, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(batchSize), static_cast<size_t>(nElectrodes)};

    std::vector<double> currentStd(batchSize * nElectrodes, 0.0);
    std::vector<size_t> currentStdShape = {static_cast<size_t>(batchSize), static_cast<size_t>(nElectrodes)};

    std::vector<double> currentSem(batchSize * nElectrodes, 0.0);
    std::vector<size_t> currentSemShape = {static_cast<size_t>(batchSize), static_cast<size_t>(nElectrodes)};
    /* 8 dimensional input space */
    std::vector<double> inputData(batchSize * nElectrodes, 0.0);
    std::vector<size_t> inputDataShape = {(size_t)batchSize, (size_t)nElectrodes};
    /* Energies  */
    /* LHC sampled voltages */
    std::vector<double> mins(nElectrodes, minVoltage);
    std::vector<double> maxs(nElectrodes, maxVoltage);
    std::vector<std::vector<double>> samples = scaledLHC(
        batchSize,
        nElectrodes,
        mins,
        maxs,
        LHCSeed 
    );
    /* Vectors for coordinates */
    std::vector<double> accCoords(batchSize*nAcceptors*2, 0.0);
    std::vector<size_t> accCoordsShape = {(size_t)batchSize, (size_t)nAcceptors, 2};
    std::vector<double> donCoords(batchSize*nDonors*2, 0.0);
    std::vector<size_t> donCoordsShape = {(size_t)batchSize, (size_t)nDonors, 2};

    /* Rate prefactors */
    std::vector<double> ratePref(batchSize*numOfSites*numOfSites, 0.0);
    std::vector<size_t> ratePrefShape = {size_t(batchSize), size_t(numOfSites), size_t(numOfSites)};

    /* Adjacency */
    std::vector<int> adjMat(batchSize * numOfSites * numOfSites, 0);
    std::vector<size_t> adjMatShape = {size_t(batchSize), size_t(numOfSites), size_t(numOfSites)};

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int _p = 0; _p < batchSize; ++_p) {
            
            uint64_t stateSeed = mix(
                uint64_t(threadBaseSeed) ^
                (uint64_t(threadID) << 32) ^
                uint64_t(_p) ^
                0xD1B54A32D192ED03ULL
            );
            uint64_t kmcSeed = mix(stateSeed ^ 0x94D049BB133111EBULL);

            /* Configuring state and simulator */
            std::vector<double> voltages = samples[_p];

            ConfigurationParams params;
            params.nAcceptors = nAcceptors;
            params.nElectrodes = nElectrodes;
            params.nDonors = nDonors;
            params.radius = radius;
            params.nu0 = nu0;
            params.a = a;
            params.T = T;
            params.energyDisorder = energyDisorder;
            params.electrodeWidth = electrodeWidth;
            params.minHopDistance = minHopDistance;
            params.maxHopDistance = maxHopDistance;
            params.Nr = Nr;
            params.Nt = Nt;
            params.distType = distType;

            uint64_t cfgSeed = mix(
                stateSeed ^ 0x94D049BB133111EBULL                
            );

            Configuration config(params, cfgSeed);

            State tempState(config, stateSeed);

            tempState.updateBoundaries(voltages);

            int s = 0;
            for (s = 0; s < numOfSites; ++s) initEn[s + _p * tempState.numOfSites] = tempState.siteEnergies[s];
            for (s = 0; s < nAcceptors; ++s) {
                randEn[s + _p * tempState.nAcceptors] = tempState.randomEnergies[s];
                accDonInt[s + _p * tempState.nAcceptors] = tempState.acceptorDonorInteraction[s];
                accAccInt[s + _p * tempState.nAcceptors] = tempState.acceptorInteraction[s];
            }

            KMCSimulator kmc(tempState, kmcSeed);
            /* Saving coords */
            int u = 0;
            for (u = 0; u < nAcceptors; ++u) {
                accCoords[(u + _p * nAcceptors)*2] = tempState.acceptorCoordinates[u * 2];
                accCoords[(u + _p * nAcceptors)*2 + 1] = tempState.acceptorCoordinates[u * 2 + 1];
            }
            for (u = 0; u < nDonors; ++u) {
                donCoords[(u + _p * nDonors)*2] = tempState.donorCoordinates[u * 2];
                donCoords[(u + _p * nDonors)*2 + 1] = tempState.donorCoordinates[u * 2 + 1];
            }

            /* Save rate prefactors into NxN shape */
            for (int k = 0; k < tempState.jaggedArrayLengths.size() - 1; ++k) {
                int start = tempState.jaggedArrayLengths[k];
                int end = tempState.jaggedArrayLengths[k+1];

                for (int l = start; l < end; ++l) {
                    ratePref[_p * numOfSites * numOfSites + k * numOfSites + tempState.neighbourIndices[l]] = kmc.constantTransitionRates[l];
                    adjMat[_p * numOfSites * numOfSites + k * numOfSites + tempState.neighbourIndices[l]] = 1;
                }
            }

            /* Equilibrate */
            kmc.simulate(tempState, eqSteps, false, false);
            tempState.resetEventCounter();
            tempState.stateTime = 0.0;

            /* Calculating output current */
            std::vector<double> averagedCurrent(nElectrodes, 0.0);
            double totalTime = 0.0;
            int intervalSteps = simSteps / numOfTasks;

            std::vector<double> meanW(nElectrodes, 0.0);
            std::vector<double> M2w(nElectrodes, 0.0);
            std::vector<double> wSum(nElectrodes, 0.0);
            std::vector<double> w2Sum(nElectrodes, 0.0);

            int intervalCount = 0;
            while (intervalCount < numOfTasks) 
            {

                double startClock = tempState.stateTime;
                kmc.simulate(tempState, intervalSteps, false, true);
                double endClock = tempState.stateTime; 

                double elapsedTime = endClock - startClock;
                totalTime += elapsedTime;
                for (int electrode_idx = 0; electrode_idx < tempState.nElectrodes; ++electrode_idx) 
                {
                    int inEvents = 0;
                    int outEvents = 0;
                    
                    for (int s = 0; s < tempState.numOfSites; ++s)
                    {
                        outEvents += tempState.eventCounter[(electrode_idx + tempState.nAcceptors) * tempState.numOfSites + s];
                        inEvents += tempState.eventCounter[tempState.numOfSites * s + (electrode_idx + tempState.nAcceptors)];                    
                    }

                    double Ii = (elapsedTime > 0) 
                        ? static_cast<double>(inEvents - outEvents) / elapsedTime : 0.0;
                    averagedCurrent[electrode_idx] += Ii * elapsedTime;

                    double w = elapsedTime;
                    double prevMean = meanW[electrode_idx];
                    wSum[electrode_idx] += w;
                    meanW[electrode_idx] += w * (Ii - meanW[electrode_idx]) / wSum[electrode_idx];
                    M2w[electrode_idx] += w * (Ii - prevMean) * (Ii - meanW[electrode_idx]);
                    w2Sum[electrode_idx] += w * w;
                }
                tempState.resetEventCounter();
                intervalCount++;
            }

            for (int electrode_idx = 0; electrode_idx < nElectrodes; ++electrode_idx)
            {
                averagedCurrent[electrode_idx] /= totalTime;
                double weightedVar = (wSum[electrode_idx] > 0.0) 
                    ? (M2w[electrode_idx] / wSum[electrode_idx]) : 0.0;
                double sampleStd = std::sqrt(weightedVar);
                double sem = (wSum[electrode_idx] > 0.0) 
                    ? (sampleStd * std::sqrt(w2Sum[electrode_idx]) / wSum[electrode_idx]) : 0.0;

                currentData[_p * nElectrodes + electrode_idx] = averagedCurrent[electrode_idx];
                currentStd[_p * nElectrodes + electrode_idx] = sampleStd;
                currentSem[_p * nElectrodes + electrode_idx] = sem;
                
            }
            for (int k = 0; k < nElectrodes; ++k) 
            {
                inputData[k + _p*nElectrodes] = voltages[k];
            }
        }
    }
    /* Input-Output */
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "w");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");
    cnpy::npz_save(file, "sampleStd", currentStd.data(), currentStdShape, "a");
    cnpy::npz_save(file, "sem", currentSem.data(), currentSemShape, "a");

    /* Coords */
    cnpy::npz_save(file, "acc_xy", accCoords.data(), accCoordsShape, "a");
    cnpy::npz_save(file, "don_xy", donCoords.data(), donCoordsShape, "a");
    
    /* Initial Energies and Single Energy contributions */
    cnpy::npz_save(file, "init_en", initEn.data(), initEnShape, "a");
    cnpy::npz_save(file, "rand_en", randEn.data(), randEnShape, "a");
    cnpy::npz_save(file, "acc_don", accDonInt.data(), accDonIntShape, "a");
    cnpy::npz_save(file, "acc_acc", accAccInt.data(), accAccIntShape, "a");

    cnpy::npz_save(file, "rate_pref", ratePref.data(), ratePrefShape, "a");
    cnpy::npz_save(file, "adj_mat", adjMat.data(), adjMatShape, "a");           
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
                if (idx == vm["outputIdx"].as<int>() || idx == vm["inputIdx"].as<int>()) {
                    std::cerr << "argParser(): Warning, c_v for input/output will be ignored";
                }
                else {
                    double v = std::stod(s.substr(eq+1));
                    voltages[idx] = v;
                }
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

    if (firstCommand == "batchFromSingleState") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("batchSize", boost::program_options::value<int>()->required())
            ("minVoltage", boost::program_options::value<double>()->required())
            ("maxVoltage", boost::program_options::value<double>()->required())
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

        batchFromSingleState(
            vm["batchSize"].as<int>(),
            vm["minVoltage"].as<double>(),
            vm["maxVoltage"].as<double>(),
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

    if (firstCommand == "batchOfMultipleStates") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("batchSize", boost::program_options::value<int>()->required())
            ("minVoltage", boost::program_options::value<double>()->required())
            ("maxVoltage", boost::program_options::value<double>()->required())
            ("nAcceptors", boost::program_options::value<int>()->required())
            ("nElectrodes", boost::program_options::value<int>()->required())
            ("nDonors", boost::program_options::value<int>()->required())
            ("radius", boost::program_options::value<double>()->required())
            ("nu0", boost::program_options::value<double>()->required())
            ("a", boost::program_options::value<double>()->required())
            ("T", boost::program_options::value<double>()->required())
            ("energyDisorder", boost::program_options::value<double>()->required())
            ("electrodeWidth", boost::program_options::value<double>()->required())
            ("minHopDistance", boost::program_options::value<double>()->required())
            ("maxHopDistance", boost::program_options::value<double>()->required())
            ("Nr", boost::program_options::value<int>()->required())
            ("Nt", boost::program_options::value<int>()->required())
            ("distType", boost::program_options::value<std::string>()->required())
            ("outputIdx", boost::program_options::value<int>()->required())
            ("eqSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simSteps", boost::program_options::value<int>()->required())
            ("numOfTasks", boost::program_options::value<int>()->default_value(100))
            ("LHCSeed", boost::program_options::value<int>()->required())
            ("threadBaseSeed", boost::program_options::value<int>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())        
            ("fileName", boost::program_options::value<std::string>()->required())
        ;
        
        boost::program_options::variables_map vm;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(options).run(),
                vm);
        boost::program_options::notify(vm);

        batchOfMultipleStates(
            vm["batchSize"].as<int>(),
            vm["minVoltage"].as<double>(),
            vm["maxVoltage"].as<double>(),
            vm["nAcceptors"].as<int>(),
            vm["nElectrodes"].as<int>(),
            vm["nDonors"].as<int>(),
            vm["radius"].as<double>(),
            vm["nu0"].as<double>(),
            vm["a"].as<double>(),
            vm["T"].as<double>(),
            vm["energyDisorder"].as<double>(),
            vm["electrodeWidth"].as<double>(),
            vm["minHopDistance"].as<double>(),
            vm["maxHopDistance"].as<double>(),
            vm["Nr"].as<int>(),
            vm["Nt"].as<int>(),
            vm["distType"].as<std::string>(),
            vm["outputIdx"].as<int>(),
            vm["eqSteps"].as<int>(),
            vm["simSteps"].as<int>(),
            vm["numOfTasks"].as<int>(),
            vm["LHCSeed"].as<int>(),
            vm["threadBaseSeed"].as<int>(),
            vm["saveFolder"].as<std::string>(),
            vm["fileName"].as<std::string>()
        );

        return 1;
    }

    return 1;
}