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

    setRandomSeed(seed);

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
    double vStep = range / static_cast<double>(numOfPoints-1);

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

void batchFromSingleState(
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
    State initState(config);
    KMCSimulator kmc(initState);
    kmc.simulate(initState, eqSteps, false, false);
    initState.resetEventCounter();
    initState.stateTime = 0.0;

    /* Vector for saving initial site energies */
    std::vector<double> initialEnergies(batchSize*initState.numOfSites, 0.0);

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int ivPoint = 0; ivPoint < batchSize; ++ivPoint) {

            int threadSeed = threadID * 100000 + threadBaseSeed + ivPoint;
            setRandomSeed(threadSeed);
            std::vector<double> voltages = samples[ivPoint];
            voltages[outputIdx] = 0.0;
            
            State equilState(initState);
            equilState.updateBoundaries(voltages);

            for (int _s = 0; _s < equilState.numOfSites; ++_s) {
                initialEnergies[_s + ivPoint*equilState.numOfSites] = equilState.siteEnergies[_s];
            }

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

            currentData[ivPoint] = averagedCurrent;
            for (int k = 0; k < config.nElectrodes; ++k) {
                inputData[k + ivPoint*config.nElectrodes] = voltages[k];
            }
        }
    }
    /* Input-Output data */
    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "w");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "a");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");

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

            ratePrefactors[initState.neighbourIndices[m] + l*initState.numOfSites] = kmc.constantTransitionRates[m];

        }
    }
    cnpy::npz_save(file, "rate_prefactors", ratePrefactors.data(), ratePrefactorsShape, "a");
}

void batchOfIndependantStates(
    int batchSize,
    double minVoltage, double maxVoltage,
    int nAcceptors, int nElectrodes, int nDonors,
    double radius, double nu0, double a, double T,
    double energyDisorder,
    double electrodeWidth,
    double minHopDistance, double maxHopDistance,
    int femRes,
    std::string distType, double epsilon,
    int inputIdx, int outputIdx,
    int eqSteps, int simSteps, int numOfTasks,
    int LHCSeed, int threadBaseSeed,
    const std::string& saveFolder, const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("[batchOfIndependantStates]: Save folder not found");
    }

    int numOfSites = nAcceptors + nElectrodes;

    std::string file = saveFolder + "/" + fileName + ".npz";
    /* Output current */
    std::vector<double> currentData(batchSize, 0.0);
    std::vector<size_t> currentDataShape = {(size_t)batchSize};
    /* 8 dimensional input space */
    std::vector<double> inputData(batchSize*nElectrodes, 0.0);
    std::vector<size_t> inputDataShape = {(size_t)batchSize, (size_t)nElectrodes};
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

    /* Vectors for saving different energy contributions */
    std::vector<double> initialEnergies(batchSize*numOfSites, 0.0);
    std::vector<double> accDonInteraction(batchSize*nAcceptors, 0.0);
    std::vector<double> accInteraction(batchSize*nAcceptors, 0.0);
    std::vector<double> randEnergies(batchSize*nAcceptors, 0.0);
    std::vector<size_t> initialEnergiesShape = {(size_t)batchSize, (size_t)numOfSites};
    std::vector<size_t> accDonInteractionShape = {(size_t)batchSize, (size_t)nAcceptors};
    std::vector<size_t> accInteractionShape = {(size_t)batchSize, (size_t)nAcceptors};
    std::vector<size_t> randEnergiesShape = {(size_t)batchSize, (size_t)nAcceptors};

    /* Vectors for coordinates */
    std::vector<double> accCoords(batchSize*nAcceptors*2, 0.0);
    std::vector<double> donCoords(batchSize*nDonors*2, 0.0);
    std::vector<size_t> accCoordsShape = {(size_t)batchSize, (size_t)nAcceptors, 2};
    std::vector<size_t> donCoordsShape = {(size_t)batchSize, (size_t)nDonors, 2};

    /* Dummy params to get electrode coordinates */
    std::vector<double> eleCoords(nElectrodes*2, 0.0);
    ConfigurationParams dummyParams;
    Configuration dummyConfig(dummyParams);
    for (int _e = 0; _e < dummyParams.electrodeData.size(); ++_e) {
        double deg = dummyConfig.electrodeData[_e].angularPosition;
        double phi = 2.0*M_PI*deg / 360.0;
        eleCoords[2*_e] =  dummyConfig.radius*std::cos(phi);                                
        eleCoords[2*_e + 1] =  dummyConfig.radius*std::sin(phi); 
    }
    std::vector<size_t> eleCoordsShape = {(size_t)nElectrodes, 2};

    /* Vector for rate prefactors */
    std::vector<double> ratePrefactors(batchSize*numOfSites*numOfSites, 0.0);
    std::vector<size_t> ratePrefactorsShape = {(size_t)batchSize, (size_t)numOfSites, (size_t)numOfSites};

    /* Neighbours */
    std::vector<int> neighbourArray(batchSize*numOfSites*numOfSites, 0);
    std::vector<size_t> neighbourArrayShape = {(size_t)batchSize, (size_t)numOfSites, (size_t)numOfSites};

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int _p = 0; _p < batchSize; ++_p) {

            int threadSeed = threadID * 100000 + threadBaseSeed + _p;
            setRandomSeed(threadSeed);
            /* Configuring state and simulator */
            std::vector<double> voltages = samples[_p];
            voltages[outputIdx] = 0.0;
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
            params.femRes = femRes;
            params.distType = distType;
            params.epsilon = epsilon;

            Configuration config(params);

            State tempState(config);

            tempState.updateBoundaries(voltages);
            KMCSimulator kmc(tempState);
            /* Saving initial energies */
            for (int _s = 0; _s < tempState.numOfSites; ++_s) {
                initialEnergies[_s + _p*tempState.numOfSites] = tempState.siteEnergies[_s];
            }
            /* Saving initial energy contributions */
            for (int _t = 0; _t < tempState.nAcceptors; ++_t) {
                accDonInteraction[_t + _p*tempState.nAcceptors] = tempState.acceptorDonorInteraction[_t];
                accInteraction[_t + _p*tempState.nAcceptors] = tempState.acceptorInteraction[_t];
                randEnergies[_t + _p*tempState.nAcceptors] = tempState.randomEnergies[_t];
            }
            /* Saving coords */
            for (int _u = 0; _u < nAcceptors; ++_u) {
                accCoords[(_u + _p*nAcceptors)*2] = tempState.acceptorCoordinates[_u*2];
                accCoords[(_u + _p*nAcceptors)*2 + 1] = tempState.acceptorCoordinates[_u*2 + 1];
            }
            for (int _u = 0; _u < nDonors; ++_u) {
                donCoords[(_u + _p*nDonors)*2] = tempState.donorCoordinates[_u*2];
                donCoords[(_u + _p*nDonors)*2 + 1] = tempState.donorCoordinates[_u*2 + 1];
            }
            /* Saving rate prefactors into NxN again */
            for (int l = 0; l < tempState.jaggedArrayLengths.size()-1; ++l) {

                int start = tempState.jaggedArrayLengths[l];
                int end = tempState.jaggedArrayLengths[l+1];

                for (int m = start; m < end; ++m) {

                    ratePrefactors[_p*numOfSites*numOfSites + l*numOfSites + tempState.neighbourIndices[m]] = kmc.constantTransitionRates[m];
                    neighbourArray[_p*numOfSites*numOfSites + l*numOfSites + tempState.neighbourIndices[m]] = 1;
                }
            }

            /* Equilibrate */
            kmc.simulate(tempState, eqSteps, false, false);
            tempState.resetEventCounter();
            tempState.stateTime = 0.0;

            /* Calculating output current */
            double averagedCurrent = 0.0;
            double totalTime = 0.0;
            int intervalSteps = simSteps / numOfTasks;
            int netEvents = 0;

            int intervalCount = 0;
            while (intervalCount < numOfTasks) {

                double startClock = tempState.stateTime;
                kmc.simulate(tempState, intervalSteps, false, true);
                double endClock = tempState.stateTime; 

                double elapsedTime = endClock - startClock;
                int inEvents = 0;
                int outEvents = 0;
                for (int i = 0; i < tempState.numOfSites; ++i) {
                    outEvents += tempState.eventCounter[(outputIdx + tempState.nAcceptors)*tempState.numOfSites + i];
                    inEvents += tempState.eventCounter[tempState.numOfSites*i + (outputIdx + tempState.nAcceptors)];
                }
                totalTime += elapsedTime;
                netEvents += inEvents-outEvents;

                tempState.resetEventCounter();

                intervalCount++;
            }

            averagedCurrent = (double)(netEvents) / totalTime;
            /* Saving output current */
            currentData[_p] = averagedCurrent;
            for (int k = 0; k < nElectrodes; ++k) {
                inputData[k + _p*nElectrodes] = voltages[k];
            }
        }
    }
    /* Input-Output data */
    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "w");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "a");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");

    /* Save coords of equilState */
    cnpy::npz_save(file, "acc_xy", accCoords.data(), accCoordsShape, "a");
    cnpy::npz_save(file, "don_xy", donCoords.data(), donCoordsShape, "a");
    cnpy::npz_save(file, "ele_xy", eleCoords.data(), eleCoordsShape, "a");

    /* Initial site energies and different energy parts */
    cnpy::npz_save(file, "init_energies", initialEnergies.data(), initialEnergiesShape, "a");
    cnpy::npz_save(file, "acc_don_int", accDonInteraction.data(), accDonInteractionShape, "a");
    cnpy::npz_save(file, "acc_acc_int", accInteraction.data(), accInteractionShape, "a");
    cnpy::npz_save(file, "rand_energies", randEnergies.data(), randEnergiesShape, "a");

    /* Constant rate prefactors */
    cnpy::npz_save(file, "rate_prefactors", ratePrefactors.data(), ratePrefactorsShape, "a");

    /* Neighbours */
    cnpy::npz_save(file, "neighbours", neighbourArray.data(), neighbourArrayShape, "a");    
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

        batchFromSingleState(
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

    if (firstCommand == "batchOfIndependantStates") {

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
            ("femRes", boost::program_options::value<int>()->required())
            ("distType", boost::program_options::value<std::string>()->required())
            ("eps", boost::program_options::value<double>()->required())
            ("inputIdx", boost::program_options::value<int>()->required())
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

        batchOfIndependantStates(
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
            vm["femRes"].as<int>(),
            vm["distType"].as<std::string>(),
            vm["eps"].as<double>(),
            vm["inputIdx"].as<int>(),
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