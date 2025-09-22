void recordDevice(
    const std::string& ID, 
    int equilibriumSteps, 
    int numOfSteps, 
    const std::string& deviceConfigs, 
    const std::string& saveFolderPath) {

    /**
     * 
     * Folder to save the batch needs to be specified.
     * 
     */

    if(saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    Simulator simulator(deviceConfigs);

    int nAcceptors = simulator.system->getAcceptorNumber();
    int nElectrodes = simulator.system->getElectrodeNumber();
    int nDonors = simulator.system->getDonorNumber();

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

    simulator.simulateNumberOfSteps(equilibriumSteps, false);
    simulator.simulateNumberOfSteps(numOfSteps, true);

    for(int i = 0; i < nAcceptors; ++i) {
        std::vector<double> coordinates = simulator.system->getAcceptorCoordinates(i);
        flattenedAcceptorCoordinates.push_back(coordinates[0]);
        flattenedAcceptorCoordinates.push_back(coordinates[1]);
    }
    for(int i = 0; i < nDonors; ++i) {
        std::vector<double> coordinates = simulator.system->getDonorCoordinates(i);
        flattenedDonorCoordinates.push_back(coordinates[0]);
        flattenedDonorCoordinates.push_back(coordinates[1]);
    }
    for(int i = 0; i < nElectrodes; ++i) {
        std::vector<double> coordinates = simulator.system->getElectrodeCoordinates(i);
        flattenedElectrodeCoordinates.push_back(coordinates[0]);
        flattenedElectrodeCoordinates.push_back(coordinates[1]);
    }

    for(int j = 0; j < nAcceptors+nElectrodes; ++j) {
        for(int i = 0; i <nAcceptors+nElectrodes; ++i) {
            flattenedEventCounts.push_back(simulator.system->getNumberOfEvents(j, i));
        }
    }

    double total_time = simulator.system->getSystemTime();

    cnpy::npz_save(deviceName, "acceptor_coordinates", flattenedAcceptorCoordinates.data(), shapeFlattenedAcceptorCoordinates, "a");
    cnpy::npz_save(deviceName, "donor_coordinates", flattenedDonorCoordinates.data(), shapeFlattenedDonorCoordinates, "a");
    cnpy::npz_save(deviceName, "electrode_coordinates", flattenedElectrodeCoordinates.data(), shapeFlattenedElectrodeCoordinates, "a");
    cnpy::npz_save(deviceName, "event_counts", flattenedEventCounts.data(), shapeFlattenedEventCounts, "a");
    cnpy::npz_save(deviceName, "device_time", &total_time, {1}, "a");
}

double IVPoint(
    std::vector<double> voltageSetting,
    int numOfDevices,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfig) {

        double averageOutputCurrent = 0.0;
        
        for (int deviceNumber = 0; deviceNumber < numOfDevices; ++deviceNumber) {
            Simulator simulator(defaultConfig);
            std::vector<double> initVoltages = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            simulator.system->updateVoltages(initVoltages);
            double current = currentFromVoltageCombination(
                simulator,
                voltageSetting,
                scanElectrodeIndex,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals,
                defaultConfig
            );

            averageOutputCurrent += current / static_cast<double>(numOfDevices);
        }

        return averageOutputCurrent;
}

double currentFromVoltageCombination(
    Simulator& simulator,
    std::vector<double> voltageSetting,
    int scanElectrodeIndex,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfig) {
        
        int nAcceptors = simulator.system->getAcceptorNumber();
        int numOfStates = simulator.system->getNumOfStates();
        simulator.system->updateVoltages(voltageSetting); 
        
        simulator.simulateNumberOfSteps(equilibriumSteps, false);

        double averageCurrent = 0.0;
        double totalTime = 0.0;
        int totalNet = 0;
        int intervalSteps = simulationSteps / numOfIntervals;
        int intervalCounter = 0;
        while(intervalCounter < numOfIntervals) {
            double startClock = simulator.system->getSystemTime();
            simulator.simulateNumberOfSteps(intervalSteps, true);
            double endClock = simulator.system->getSystemTime();

            double elapsedTime = endClock - startClock;
            int inCounts = 0;
            int outCounts = 0;
            for (int i = 0; i < numOfStates; ++i) {
                outCounts += simulator.system->getNumberOfEvents(nAcceptors+scanElectrodeIndex, i);
                inCounts += simulator.system->getNumberOfEvents(i, nAcceptors+scanElectrodeIndex); 
            }
            totalTime += elapsedTime;
            totalNet += inCounts-outCounts;
            
            simulator.system->resetEventCounts();
            intervalCounter++;
        }

        averageCurrent = static_cast<double>(totalNet) / totalTime;

        return averageCurrent;
}

void createBatchOfSingleSystem(
    int batchSize, 
    int outputElectrodeIndex,
    double minVoltage, 
    double maxVoltage,
    int equilibriumSteps,
    int simulationSteps,
    int numOfIntervals,
    const std::string& defaultConfigs, 
    const std::string& saveFolderPath, 
    const std::string& batchID) {

    if (saveFolderPath.empty()) {
        throw std::invalid_argument("No save folder specified !");
    }

    std::string fileName = saveFolderPath + "/batch_" + batchID + ".npz";

    std::vector<int> systemElectrodes = {0, 1, 2, 3, 4, 5, 6, 7};
    int numOfElectrodes = systemElectrodes.size();

    std::vector<double> inputs(batchSize*numOfElectrodes, 0.0);
    std::vector<double> outputs(batchSize, 0.0);   
    std::vector<size_t> shapeInputs = {static_cast<size_t>(batchSize), static_cast<size_t>(numOfElectrodes)}; 
    std::vector<size_t> shapeOutputs = {static_cast<size_t>(batchSize)};  
    
    double temperature = 0.0;
    double sigma = 0.0;
    int nDonors = 0;

    #pragma omp parallel 
    {   
        std::mt19937 rng(std::random_device{}() + omp_get_thread_num());
        std::uniform_real_distribution<double> uni(minVoltage, maxVoltage);

        Simulator simulator(defaultConfigs);
        
        #pragma omp for schedule(static)
        for (int batch = 0; batch < batchSize; ++batch) {
            std::vector<double> voltages = sampleVoltageSetting(numOfElectrodes, -1.5, 1.5);//{0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0};
            
            if (batch == 0) {
                temperature = simulator.system->T;
                sigma = simulator.system->energyDisorder;
                nDonors = simulator.system->nDonors;
            }
            std::vector<double> initVoltages = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            simulator.system->updateVoltages(initVoltages);

            int nAcceptors = simulator.system->getAcceptorNumber();
            int numOfStates = simulator.system->getNumOfStates();
            
            voltages[(outputElectrodeIndex+1) % 8] = uni(rng);
            voltages[outputElectrodeIndex] = 0.0;

            for (int i = 0; i < numOfElectrodes; ++i) {
                inputs[batch*numOfElectrodes + i] = voltages[i];
            }

            double averageOutputCurrent = currentFromVoltageCombination(
                simulator,
                voltages,
                outputElectrodeIndex,
                equilibriumSteps,
                simulationSteps,
                numOfIntervals, 
                defaultConfigs
            );
            
            outputs[batch] = averageOutputCurrent;
        }
    }

    cnpy::npz_save(fileName, "ID", &batchID, {1}, "w");
    cnpy::npz_save(fileName, "inputs", inputs.data(), shapeInputs, "a");
    cnpy::npz_save(fileName, "outputs", outputs.data(), shapeOutputs, "a");
    cnpy::npz_save(fileName, "T/K", &temperature, {1}, "a");
    cnpy::npz_save(fileName, "sigma/kbT", &sigma, {1}, "a");
    cnpy::npz_save(fileName, "nDonors", &nDonors, {1}, "a");
}

int argParser(int argc, char* argv[]) {

    boost::program_options::options_description globalOptions(" ");

    globalOptions.add_options()
        ("help, h", "help message")
        ("command", boost::program_options::value<std::string>(), "command to run")
    ;
    
    boost::program_options::positional_options_description position;
    position.add("command", 1);

    auto parsedCommand = boost::program_options::command_line_parser(argc, argv).options(globalOptions)
                                                                                .positional(position)
                                                                                .allow_unregistered()
                                                                                .run();
    
    boost::program_options::variables_map vm;
    boost::program_options::store(parsedCommand, vm);

    if (vm.count("help") || !vm.count("command")) {
        std::cout << globalOptions << "\n"
                  << "Allowed commands:\n"
                  << " singleRun --equilibriumSteps <int> --simulationSteps <int> --deviceName <string>\n"
                  << " batchRun --batchSize <int> --equilibriumSteps <int> --simulationSteps <int> --batchName <string>\n";
        return 0;
    }

    std::string firstCommand = vm["command"].as<std::string>();
    std::vector<std::string> remainingCommand = 
        boost::program_options::collect_unrecognized(
            parsedCommand.options, 
            boost::program_options::include_positional);

    remainingCommand.erase(remainingCommand.begin());

    if (firstCommand == "singleRun") {
        boost::program_options::options_description singleRunOptions("Single run options");
        singleRunOptions.add_options()
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("deviceName", boost::program_options::value<std::string>()->required())
        ;
        
        boost::program_options::variables_map singleRunVM;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(singleRunOptions).run(),
                singleRunVM);
        boost::program_options::notify(singleRunVM);

        recordDevice(
            singleRunVM["deviceName"].as<std::string>(),
            singleRunVM["equilibriumSteps"].as<int>(),
            singleRunVM["simulationSteps"].as<int>(),
            "default_configs",
            "currentData"
        );
        
        return 1;
    }

    else if (firstCommand == "batchRun") {
        boost::program_options::options_description batchRunOptions("Batch run options");
        batchRunOptions.add_options()
            ("batchSize", boost::program_options::value<int>()->default_value(512))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("batchName", boost::program_options::value<std::string>()->required())
            ("saveFolderPath", boost::program_options::value<std::string>()->required())
            ("min", boost::program_options::value<double>()->default_value(-1.5))
            ("max", boost::program_options::value<double>()->default_value(1.5))
        ;
        
        boost::program_options::variables_map batchRunVM;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(batchRunOptions).run(),
                batchRunVM);
        boost::program_options::notify(batchRunVM);
        
        createBatchOfSingleSystem(
            batchRunVM["batchSize"].as<int>(),
            0,
            -1.5,
            1.5,
            batchRunVM["equilibriumSteps"].as<int>(),
            batchRunVM["simulationSteps"].as<int>(),
            10, 
            "default_configs",
            batchRunVM["saveFolderPath"].as<std::string>(),
            batchRunVM["batchName"].as<std::string>()
        );

        return 1;
    }

    else if (firstCommand == "voltageSweep") {
        boost::program_options::options_description sweepOptions("Sweep run options");
        sweepOptions.add_options()
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("fileName", boost::program_options::value<std::string>()->required())
            ("min", boost::program_options::value<double>()->default_value(-1.5))
            ("max", boost::program_options::value<double>()->default_value(1.5))
        ;

        boost::program_options::variables_map sweepVM;
        boost::program_options::store(
            boost::program_options::command_line_parser(
                remainingCommand).options(sweepOptions).run(),
                sweepVM);
        boost::program_options::notify(sweepVM);

        std::string saveFolderPath = "currentData/sweep_" + sweepVM["fileName"].as<std::string>() + ".npz";
        
        int outputElectrodeIndex = 0;
        int inputElectrodeIndex = 5;
        std::vector<double> inputs(sweepVM["numOfPoints"].as<int>(), 0.0);
        std::vector<double> outputs(sweepVM["numOfPoints"].as<int>(), 0.0);
        
        std::vector<double> controlVoltages = sampleVoltageSetting(8, sweepVM["min"].as<double>(), sweepVM["max"].as<double>());
        controlVoltages[outputElectrodeIndex] = 0.0;

        std::vector<double> initSetting {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        Simulator sim("default_configs");
        sim.system->updateVoltages(initSetting);   
        double averageCurrent = 0.0;

        for (int i = 0; i < sweepVM["numOfPoints"].as<int>(); ++i) {
            controlVoltages[inputElectrodeIndex] = sampleFromUniformDistribution(sweepVM["min"].as<double>(), sweepVM["max"].as<double>());
            sim.system->updateVoltages(controlVoltages);
            averageCurrent = currentFromVoltageCombination(
                sim,
                controlVoltages,
                outputElectrodeIndex,
                sweepVM["equilibriumSteps"].as<int>(),
                sweepVM["simulationSteps"].as<int>(),
                1,
                "default_configs"
            );

            inputs[i] = controlVoltages[inputElectrodeIndex];
            outputs[i] = averageCurrent;
        }

        controlVoltages[outputElectrodeIndex] = 0.0;
        controlVoltages[inputElectrodeIndex] = 0.0;

        std::vector<size_t> controlVoltagesShape = {static_cast<size_t>(controlVoltages.size())};
        std::vector<size_t> dataShape = {static_cast<size_t>(sweepVM["numOfPoints"].as<int>())};
        
        cnpy::npz_save(saveFolderPath, "input_idx", &inputElectrodeIndex, {1}, "w");
        cnpy::npz_save(saveFolderPath, "output_idx", &outputElectrodeIndex, {1}, "a");
        cnpy::npz_save(saveFolderPath, "control_voltages", controlVoltages.data(), controlVoltagesShape, "a");
        cnpy::npz_save(saveFolderPath, "inputs", inputs.data(), dataShape, "a");
        cnpy::npz_save(saveFolderPath, "outputs", outputs.data(), dataShape, "a");
    }
    
    else {
        std::cerr << "Unknown command: " << firstCommand << "\n";
        return -1;
    }

    return 0;
}

State::State()
    : nAcceptors(200)
    , nDonors(3)
    , nElectrodes(8)
    , numOfSites(nAcceptors + nElectrodes)
    , radius(150.0)
    , nu0(1.0)
    , a(20.0)
    , T(77.0)
    , kbT(kb*T)
    , energyDisorder(0.05*e / kbT)
    , minHopDistance(3.0)
    , maxHopDistance(60.0)
    , noDimension(true)
{
    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    currentOccupation.resize(nAcceptors, 0);
    initialOccupation.resize(nAcceptors, 0);
    randomEnergies.resize(nAcceptors, 0.0);
    acceptorDonorInteraction.resize(nAcceptors, 0.0);
    acceptorInteraction.resize(nAcceptors*nAcceptors, 0.0);
    initialSiteEnergies.resize(nAcceptors+nElectrodes, 0.0);
    initialPotential.resize(nAcceptors+nElectrodes, 0.0);
    currentPotential.resize(nAcceptors+nElectrodes, 0.0);
    siteEnergies.resize(nAcceptors+nElectrodes, 0.0);
    eventCounter.resize(numOfSites*numOfSites, 0);

    initRandomState();
}

void State::initRandomState() {

    R = std::sqrt(M_PI*radius*radius / static_cast<double>(nAcceptors));

    if (noDimension) {
        radius = radius / R;
    }

    for (int i = 0; i < nAcceptors; ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        acceptorCoordinates[i*2] = randomR*std::cos(randomPhi);
        acceptorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }

    for (int i = 0; i < nDonors;  ++i) {
        double randomPhi = 2.0*M_PI*sampleFromUniformDistribution(0.0, 1.0);
        double randomR = radius*std::sqrt(sampleFromUniformDistribution(0.0, 1.0));
        donorCoordinates[i*2] = randomR*std::cos(randomPhi);
        donorCoordinates[i*2 + 1] = randomR*std::sin(randomPhi);
    }  

    electrodeData.resize(8);

    std::vector<double> defaultPositions = {
        0.0,
        45.,
        90.,
        135.,
        180.,
        225.,
        270.,
        315.,
        360.
    };

    for (int i = 0; i < electrodeData.size(); ++i) {
        electrodeData[i].angularPosition = defaultPositions[i];
        electrodeData[i].voltage = -1.5 + 3.0*randomDouble01();
    }
}

void KMCSimulator::initKMCSimulator(State& state) {

    numOfNeighbors = state.numOfNeighbours;
    jaggedArrayLengths = state.jaggedArrayLengths;
    neighborIndices = state.neighbourIndices;
    totalNumOfEvents = state.totalNumOfEvents;
    lastHopIndices.resize(2, 0);

    constantTransitionRates.resize(2*totalNumOfEvents, 0.0);
    dynamicalTransitionRates.resize(2*totalNumOfEvents, 0.0);
    aggregatedTransitionRates.resize(2*totalNumOfEvents, 0.0);

    std::vector<int> writePtr(state.numOfSites);
    for (int i = 0; i < state.numOfSites; ++i) {
        writePtr[i] = jaggedArrayLengths[i];
    }

    for (int i = 0; i < state.numOfSites; ++i) {
        for (int j = i+1; j < state.numOfSites; ++j) {
        double distance =  state.distanceMatrix[i*state.numOfSites + j];
            if (distance > state.minHopDistance && distance < state.maxHopDistance) {
                int indexIJ = writePtr[i]++;
                int indexJI = writePtr[j]++;
                constantTransitionRates[indexIJ] = state.nu0*fastExp(-2.0*distance / state.a);
                constantTransitionRates[indexJI] = state.nu0*fastExp(-2.0*distance / state.a);
            }
        }
    }
}

void State::initStateFromConfig(Configuration& config) {

    nAcceptors = config.nAcceptors;
    nDonors = config.nDonors;
    nElectrodes = config.nElectrodes;
    numOfSites = config.numOfSites;
    radius = config.radius;
    nu0 = config.nu0;
    a = config.a;
    T = config.T;
    energyDisorder = config.energyDisorder;
    R = config.R;
    A0 = config.A0;
    electrodeWidth = config.electrodeWidth;
    minHopDistance = config.minHopDistance;
    maxHopDistance = config.maxHopDistance;

    acceptorCoordinates.resize(2*nAcceptors, 0.0);
    donorCoordinates.resize(2*nDonors, 0.0);
    electrodeCoordinates.resize(2*nElectrodes, 0.0);
    
    distanceMatrix.resize(numOfSites*numOfSites, 0.0);
    inverseAcceptorDistances.resize(nAcceptors*nAcceptors, 0.0);
    currentOccupation.resize(nAcceptors, 0);
    initialOccupation.resize(nAcceptors, 0);
    randomEnergies.resize(nAcceptors, 0.0);
    acceptorDonorInteraction.resize(nAcceptors, 0.0);
    acceptorInteraction.resize(nAcceptors*nAcceptors, 0.0);
    initialSiteEnergies.resize(nAcceptors+nElectrodes, 0.0);
    initialPotential.resize(nAcceptors+nElectrodes, 0.0);
    currentPotential.resize(nAcceptors+nElectrodes, 0.0);
    siteEnergies.resize(nAcceptors+nElectrodes, 0.0);   
    eventCounter.resize(numOfSites*numOfSites, 0); 
}

// PARALLEL LOOP EXAMPLE

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

            //std::cout << "Finished batch#" << _batch << "\n";
            //std::cout << minVoltage + _batch*vStep << "\n";
        }
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
    std::cout << "Entering lineSweep()" << "\n";
    if (saveFolderPath.empty()) {
        throw std::invalid_argument("lineSweep(): No save folder found");
    }
    std::string filePath = saveFolderPath + "/" + fileName + ".npz";
    std::cout << "N: " << N << "\n"; 
    std::cout << "sampleSize: " << sampleSize << "\n";
    std::cout << "numOfPoints: " << numOfPoints << "\n";
    std::cout << "OutputLength: " << N*numOfPoints << "\n";
    int femRes = 1e5;
    int numOfElectrodes = 8;

    double range = vMax - vMin;
    double vStep = range / static_cast<double>(numOfPoints - 1);

    double paramStep = (varParamMax - varParamMin) / static_cast<double>(N-1);

    std::vector<double> outputs(N*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(N), static_cast<size_t>(numOfPoints)};
    std::cout << "Entering parallel loop:" << "\n";
    #pragma omp parallel
    {   
        #pragma omp for
        for (int i = 0; i < N; ++i) {
            std::cout << "i: " << i << "\n";
            for (int s = 0; s < sampleSize; ++s) {
                std::cout << "s: " << s << "\n";
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
                    std::cout << "_b: " << _b << "\n";
                    newBoundaries[_b] = voltages[_b];
                }
                newBoundaries[outputIdx] = 0.6;

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
                    std::cout << "_v: " << _v << "\n";
                    outputs[i*numOfPoints + _v] += current / static_cast<double>(sampleSize);
                }
            }
        }
    }
    cnpy::npz_save(filePath, "ID", &fileName, {1}, "w");
    cnpy::npz_save(filePath, "out", outputs.data(), outputShape, "a");
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

        std::cout << "Entering argParser with lineSweep arg" << "\n";
        
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
        std::cout << "Voltages in argParser" << "\n";
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
        std::cout << "Calling lineSweep from argParser:" << "\n";
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
);

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
);

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
);

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
);

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
);

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

void singleIVCurve (
    int numOfPoints,
    int inputIdx,
    int outputIdx,
    double minVoltage,
    double maxVoltage,
    int eqSteps,
    int simSteps,
    int numIntervals,
    int seed,
    std::vector<double> controlVoltages,
    const std::string& cfg,
    const std::string& acceptorCfg,
    const std::string& donorCfg,
    const std::string& electrodeCfg,
    const std::string& saveFolder,
    const std::string& fileName
) {

    if (saveFolder.empty()) {
        throw std::invalid_argument("singleIVCurve(): Save folder not found");
    }
    std::string file = saveFolder + "/" + fileName + ".npz";

    std::vector<double> currentData(numOfSamples*numOfPoints, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(numOfSamples), static_cast<size_t>(numOfPoints)};

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
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        int threadSeed = seed + threadID * 1000000;

        #pragma omp for
        for (int s = 0; s < numOfSamples; ++s) {

            int sampleSeed = threadSeed + s;
            setRandomSeed(sampleSeed);

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

            std::vector<double> newBoundaries(controlVoltages.size(), 0.0);
            newBoundaries[outputIdx] = 0.0;
            for (int v = 0; v < numOfPoints; v++) {

                newBoundaries[inputIdx] = minVoltage + v*vStep;
                state.updateBoundaries(newBoundaries, fem);

                double averagedCurrent = calculateCurrent(
                    state,
                    kmc,
                    outputIdx,
                    eqSteps,
                    simSteps,
                    numIntervals
                );

                currentData[s*numOfPoints + v] = averagedCurrent;
            }
        }
    }

    cnpy::npz_save(file, "current", currentData.data(), currentDataShape, "w");
    cnpy::npz_save(file, "control", controlData.data(), controlDataShape, "a");
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
    const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("singleStatebatch: Save folder not found");
    }
    std::string file = saveFolder + "/" + fileName + ".npz";

    int femResolution = 1e5;
    int numOfElectrodes = 8;

    std::vector<double> outputs(batchSize*numOfPoints, 0.0);
    std::vector<size_t> outputShape = {static_cast<size_t>(batchSize), static_cast<size_t>(numOfPoints)}; 

    double range = maxVoltage - minVoltage;
    double vStep = range / (numOfPoints - 1);

    //std::vector<double> mins = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
    //std::vector<double> maxs = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    std::vector<double> mins(8, 0.0);
    std::vector<double> maxs(8, 0.0);
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
    
    cnpy::npz_save(file, "out", outputs.data(), outputShape, "w");
}


    if (firstCommand == "batchRun") {

        boost::program_options::options_description options("Batch run options");
        options.add_options()
            ("batchSize", boost::program_options::value<int>()->required())
            ("numOfPoints", boost::program_options::value<int>()->default_value(100))
            ("numOfCurves", boost::program_options::value<int>()->required())
            ("inputIdx", boost::program_options::value<int>()->required())
            ("outputIdx", boost::program_options::value<int>()->required())
            ("minVoltage", boost::program_options::value<double>()->required())
            ("maxVoltage", boost::program_options::value<double>()->required())
            ("equilibriumSteps", boost::program_options::value<int>()->default_value(1e4))
            ("simulationSteps", boost::program_options::value<int>()->required())
            ("numOfIntervals", boost::program_options::value<int>()->default_value(100))
            ("seed", boost::program_options::value<int>()->required())
            ("cfg", boost::program_options::value<std::string>()->required())
            ("acceptorCfg", boost::program_options::value<std::string>()->required())
            ("donorCfg", boost::program_options::value<std::string>()->required())
            ("electrodeCfg", boost::program_options::value<std::string>()->required())
            ("saveFolder", boost::program_options::value<std::string>()->required())        
            ("fileName", boost::program_options::value<std::string>()->required())
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
            vm["minVoltage"].as<double>(),
            vm["maxVoltage"].as<double>(),
            vm["equilibriumSteps"].as<int>(),
            vm["simulationSteps"].as<int>(),
            vm["numOfIntervals"].as<int>(),
            vm["seed"].as<int>(),
            vm["cfg"].as<std::string>(),
            vm["acceptorCfg"].as<std::string>(),
            vm["donorCfg"].as<std::string>(),
            vm["electrodeCfg"].as<std::string>(),
            vm["saveFolder"].as<std::string>(),
            vm["fileName"].as<std::string>()
        );

        return 1;
    }

    if (firstCommand == "singleCurve") {

        boost::program_options::options_description options("Find control voltages options");
        options.add_options()
            ("numOfPoints", boost::program_options::value<int>()->required())
            ("numOfSamples", boost::program_options::value<int>()->required())
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
        std::cout << vm["cfg"].as<std::string>() << "\n";
        singleIVCurve(
            vm["numOfPoints"].as<int>(),
            vm["numOfSamples"].as<int>(),
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

    State& State::operator=(State const& other) {
    if (this == &other) return *this;

    acceptorCoordinates = other.acceptorCoordinates;
    donorCoordinates = other.donorCoordinates;
    electrodeCoordinates = other.electrodeCoordinates;
    nAcceptors = other.nAcceptors;
    nDonors = other.nDonors;
    nElectrodes = other.nElectrodes;
    numOfSites = other.numOfSites;
    radius = other.radius;
    nu0 = other.nu0;
    a = other.a;
    T = other.T;
    kbT = other.kbT;
    energyDisorder = other.energyDisorder;
    R = other.R;
    A0 = other.A0;
    electrodeWidth = other.electrodeWidth;
    minHopDistance = other.minHopDistance;
    maxHopDistance = other.maxHopDistance;
    electrodeData = other.electrodeData;

    distanceMatrix = other.distanceMatrix;
    inverseAcceptorDistances = other.inverseAcceptorDistances;
    currentOccupation = other.currentOccupation;
    initialOccupation = other.initialOccupation;
    randomEnergies = other.randomEnergies;
    acceptorDonorInteraction = other.acceptorDonorInteraction;
    acceptorInteraction = other.acceptorInteraction;
    initialSiteEnergies = other.initialSiteEnergies;
    initialPotential = other.initialPotential;
    currentPotential = other.currentPotential;
    siteEnergies = other.siteEnergies;
    eventCounter = other.eventCounter;
    jaggedArrayLengths = other.jaggedArrayLengths;
    neighbourIndices = other.neighbourIndices;
    numOfNeighbours = other.numOfNeighbours;
    totalNumOfEvents = other.totalNumOfEvents;
    stateTime = other.stateTime;
    femRes = other.femRes;

    return *this;
}

while (true) {
                        double x = normalDist(0.0, 1.0);
                        double y = normalDist(0.0, 1.0);

                        if (std::sqrt(x*x + y*y) <= radius) {
                            double randomPhi = 2.0*M_PI*randomDouble01();
                            acceptorCoords.push_back(radius*std::cos(randomPhi));
                            acceptorCoords.push_back(radius*std::sin(randomPhi));
                            break;
                        }
                    }


else {

        if (type == "uniform") {

            auto now = std::chrono::high_resolution_clock::now();
            auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
            setRandomSeed(static_cast<long int>(now_ns));

            for (int i = 0; i < nAcceptors; ++i) {
                double randomPhi = 2.0*M_PI*randomDouble01();
                double randomR = radius*std::sqrt(randomDouble01());
                acceptorCoords.push_back(randomR*std::cos(randomPhi));
                acceptorCoords.push_back(randomR*std::sin(randomPhi));
            }

            for (int i = 0; i < nDonors; ++i) {
                double randomPhi = 2.0*M_PI*randomDouble01();
                double randomR = radius*std::sqrt(randomDouble01());
                donorCoords.push_back(randomR*std::cos(randomPhi));
                donorCoords.push_back(randomR*std::sin(randomPhi));
            }
        }

        else if (type == "mixed") {

            if (epsilon < 0.0 || epsilon > 1.0) {
                throw std::invalid_argument("Configuration: epsilon must be in [0,1]");
            }

            auto now = std::chrono::high_resolution_clock::now();
            auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
            setRandomSeed(static_cast<long int>(now_ns)); 

            for (int i = 0; i < nAcceptors; ++i) {
                
                double u01 = randomDouble01();

                if (u01 < (1 - epsilon)) {
                    double randomPhi = 2.0*M_PI*randomDouble01();
                    double randomR = radius*std::sqrt(randomDouble01());
                    acceptorCoords.push_back(randomR*std::cos(randomPhi));
                    acceptorCoords.push_back(randomR*std::sin(randomPhi));
                }
                else {
                    double stdScaled = radius / 2.5;
                    std::vector<double> coords = sample_truncated_gaussian_reject(stdScaled, radius);

                    double _r = std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
                    double randomPhi = 2.0*M_PI*randomDouble01();

                    acceptorCoords.push_back(coords[0]);
                    acceptorCoords.push_back(coords[1]);
                }
            }

            for (int i = 0; i < nDonors; ++i) {
                double randomPhi = 2.0*M_PI*randomDouble01();
                double randomR = radius*std::sqrt(randomDouble01());
                donorCoords.push_back(randomR*std::cos(randomPhi));
                donorCoords.push_back(randomR*std::sin(randomPhi));
            }
        }
    }

    void batchOfIVPointsWithDistParam(
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
    const std::string& distType,
    const std::string& saveFolder,
    const std::string& fileName
) {
    if (saveFolder.empty()) {
        throw std::invalid_argument("[batchOfIVPointsWithDistParam]: Save folder not found");
    }

    ConfigurationParams cfgParams;

    std::string file = saveFolder + "/" + fileName + ".npz";

    std::vector<double> currentData(batchSize, 0.0);
    std::vector<size_t> currentDataShape = {static_cast<size_t>(batchSize)};
    /* 8 electrode params + 1 (eps) */
    int numParams = 9;

    std::vector<double> inputData(batchSize*numParams, 0.0);
    std::vector<size_t> inputDataShape = {static_cast<size_t>(batchSize), static_cast<size_t>(numParams)};

    /*coordinate data*/
    std::vector<double> accCoords(batchSize*nAcceptors*2);
    std::vector<size_t> accCoordsShape = {static_cast<size_t>(batchSize), static_cast<size_t>(nAcceptors), 2};

    std::vector<double> donCoords(batchSize*nDonors*2);
    std::vector<size_t> donCoordsShape = {static_cast<size_t>(batchSize), static_cast<size_t>(nDonors), 2};

    /*boundary vls are now for voltages + eps*/
    std::vector<double> minParamValues(numParams, -1.5);
    std::vector<double> maxParamValues(numParams, 1.5);
    /*sets boundary vls for eps*/
    minParamValues[8] = 0.0;
    maxParamValues[8] = 1.0;

    std::vector<std::vector<double>> paramSamples = scaledLHC(
        batchSize,
        numParams,
        minParamValues,
        maxParamValues,
        LHCSeed
    );

    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        #pragma omp for
        for (int _p = 0; _p < batchSize; ++_p) {

            int threadSeed = threadID * 100000 + threadBaseSeed + _p;
            setRandomSeed(threadSeed);

            std::vector<double> params = paramSamples[_p];
            double epsParam = params[numParams-1];

            Configuration config(
                cfgParams
            );

            State equilState(config);
            KMCSimulator kmc(equilState);
            kmc.simulate(equilState, eqSteps, false, false);

            equilState.resetEventCounter();
            equilState.stateTime = 0.0;

            std::vector<double> voltages = std::vector<double>(params.begin(), params.end() - 1);

            voltages[outputIdx] = 0.0;

            double averagedCurrent = singleIVPoint(
                equilState,
                outputIdx,
                numOfTasks,
                simSteps,
                voltages
            );

            currentData[_p] = averagedCurrent;
            for (int k = 0; k < numParams; ++k) {
                inputData[k + _p*numParams] = params[k];
            }     
            
            for (int l = 0; l < equilState.nAcceptors; ++l) {
                accCoords[2*l + _p*equilState.nAcceptors*2] = equilState.acceptorCoordinates[l*2];
                accCoords[2*l + _p*equilState.nAcceptors*2 + 1] = equilState.acceptorCoordinates[l*2 + 1];
            }

            for (int m = 0; m < equilState.nDonors; ++m) {
                donCoords[2*m + _p*equilState.nDonors*2] = equilState.donorCoordinates[m*2];
                donCoords[2*m + _p*equilState.nDonors*2 + 1] = equilState.donorCoordinates[m*2 + 1];
            }
        }
    }
    /*input-output data*/
    cnpy::npz_save(file, "inputIdx", &inputIdx, {1}, "w");
    cnpy::npz_save(file, "outputIdx", &outputIdx, {1}, "a");
    cnpy::npz_save(file, "currents", currentData.data(), currentDataShape, "a");
    cnpy::npz_save(file, "inputs", inputData.data(), inputDataShape, "a");
    /*coordinate data*/
    cnpy::npz_save(file, "acc_xy", accCoords.data(), accCoordsShape, "a");
    cnpy::npz_save(file, "don_xy", donCoords.data(), donCoordsShape, "a");
}

if (firstCommand == "batch_with_dist_param") {

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
            ("distType", boost::program_options::value<std::string>()->required())
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

        batchOfIVPointsWithDistParam(
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
            vm["distType"].as<std::string>(),
            vm["cfg"].as<std::string>(),
            vm["accCfg"].as<std::string>(),
            vm["donCfg"].as<std::string>(),
            vm["eleCfg"].as<std::string>(),
            vm["saveFolder"].as<std::string>(),
            vm["fileName"].as<std::string>()
        );

        return 1;
    }

    void batchOfIVPointsWithDistParam(
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
    const std::string& distType,
    const std::string& cfg,
    const std::string& accCfg,
    const std::string& donCfg,
    const std::string& eleCfg,
    const std::string& saveFolder,
    const std::string& fileName
);

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

    double randomUniform01() {

    static thread_local std::mt19937 generator;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(generator);
}

double randomNormal(double mean, double stdDev) {

    static thread_local std::mt19937 generator;
    std::normal_distribution<double> dist(mean, stdDev);

    return dist(generator);
}

double randomUniform01();

double randomNormal(double mean, double stdDev); 

/* else if (distType == "mixed") {

        if (epsilon < 0.0 || epsilon > 1.0) {
            throw std::invalid_argument("[Configuration]: invalid epsilon value");
        }

        auto now = std::chrono::high_resolution_clock::now();
        auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        setRandomSeed(static_cast<long int>(now_ns)); 

        for (int i = 0; i < nAcceptors; ++i) {
            
            double u01 = randomDouble01();

            if (u01 < (1 - epsilon)) {
                double randomPhi = 2.0*M_PI*randomDouble01();
                double randomR = radius*std::sqrt(randomDouble01());
                acceptorCoords.push_back(randomR*std::cos(randomPhi));
                acceptorCoords.push_back(randomR*std::sin(randomPhi));
            }
            else {
                double stdScaled = radius / 2.5;
                std::vector<double> coords = sample_truncated_gaussian_reject(stdScaled, radius);

                double _r = std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
                double randomPhi = 2.0*M_PI*randomDouble01();

                acceptorCoords.push_back(coords[0]);
                acceptorCoords.push_back(coords[1]);
            }
        }

        for (int i = 0; i < nDonors; ++i) {
            double randomPhi = 2.0*M_PI*randomDouble01();
            double randomR = radius*std::sqrt(randomDouble01());
            donorCoords.push_back(randomR*std::cos(randomPhi));
            donorCoords.push_back(randomR*std::sin(randomPhi));
        }
    } */

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

/* int num_electrodes = 8;
    int batch_size = 200;
    std::vector<double> mins(num_electrodes, -1.5);
    std::vector<double> maxs(num_electrodes, 1.5);
    std::vector<std::vector<double>> samples = scaledLHC(
        batch_size,
        num_electrodes,
        mins,
        maxs, 
        32
    );
    double radius = 150.0;
    double electrodeWidth = 60.0;
    int fem_res = 100000;
    for (int i = 0; i < batch_size; ++i) {
        
        std::unique_ptr<FiniteElementeCircle> femSolver;
    
        auto femCircle = std::make_unique<FiniteElementeCircle>(
            radius, fem_res
        );
        ConfigurationParams cfgParams;
        Configuration cfg(cfgParams, 32);

        auto electrodeData = cfgParams.electrodeData;

        for (int k = 0; k < electrodeData.size(); ++k) {
            femCircle->setElectrode(
                0.0,
                electrodeData[i].angularPosition/360.0 * 2.0*M_PI - 0.5*electrodeWidth / radius,
                electrodeData[i].angularPosition/360.0 * 2.0*M_PI + 0.5*electrodeWidth / radius
            );
        }

        femSolver = std::move(femCircle);
        femSolver->initRun();
    }

    std::cout << "Hello" << "\n"; */

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


    class LiteLaplaceCircle {
public:
    // Construct a polar grid on a disk of radius R.
    // Nr: number of radial nodes (>=2), Ntheta: number of angular nodes (>=8).
    // useSOR: if true, Successive Over-Relaxation with omega (1<omega<2).
    explicit LiteLaplaceCircle(double R, int Nr, int Ntheta,
                               bool useSOR = true, double omega = 1.8)
        : R_(R), Nr_(Nr), Nt_(Ntheta),
          useSOR_(useSOR), omega_(omega)
    {
        if (R_ <= 0)               throw std::invalid_argument("R must be > 0");
        if (Nr_ < 2)               throw std::invalid_argument("Nr must be >= 2");
        if (Nt_ < 8)               throw std::invalid_argument("Ntheta must be >= 8");
        if (useSOR_ && (omega_ <= 1.0 || omega_ >= 2.0))
            throw std::invalid_argument("omega in (1,2) for SOR");

        dr_ = R_ / (Nr_ - 1);
        dth_ = 2.0 * M_PI / Nt_;

        phi_.assign(Nr_ * Nt_, 0.0);
        phi_new_.assign(Nr_ * Nt_, 0.0);

        // Dirichlet mask only lives on the outer ring (r=R)
        outer_dirichlet_mask_.assign(Nt_, false);
        outer_dirichlet_value_.assign(Nt_, 0.0);
    }

    // Resets the field to zeros (handy between tests)
    void clearSolution(double value = 0.0) {
        std::fill(phi_.begin(), phi_.end(), value);
    }

    // Define an electrode arc on the outer boundary with a constant voltage.
    // begin/end are in radians; any order is allowed and wrapped to [0,2pi).
    void setElectrode(double voltage, double begin, double end) {
        auto wrap = [](double a) {
            double t = std::fmod(a, 2.0*M_PI);
            return (t < 0.0) ? t + 2.0*M_PI : t;
        };
        begin = wrap(begin);
        end   = wrap(end);

        auto mark_index = [&](int j) {
            outer_dirichlet_mask_[j]  = true;
            outer_dirichlet_value_[j] = voltage;
        };

        if (begin <= end) {
            int j0 = thetaToIndex(begin);
            int j1 = thetaToIndex(end);
            // cover inclusive range on the ring
            for (int j = j0; j <= j1; ++j) mark_index(j % Nt_);
        } else {
            // wrapped arc
            int j0 = thetaToIndex(begin);
            int j1 = thetaToIndex(end);
            for (int j = j0; j < j0 + Nt_; ++j) {
                int jj = j % Nt_;
                mark_index(jj);
                if (jj == j1) break;
            }
        }
        applyOuterBoundaryToSolution();
        electrodes_.push_back(snapshotElectrodeIndices());
    }

    // Change voltage of an existing electrode (by creation order).
    void updateElectrodeVoltage(int electrodeIndex, double voltage) {
        if (electrodeIndex < 0 || electrodeIndex >= (int)electrodes_.size())
            throw std::out_of_range("Bad electrode index");
        // Clear previous values for this electrode, then re-apply at new voltage.
        for (int j : electrodes_[electrodeIndex]) {
            outer_dirichlet_mask_[j]  = true;
            outer_dirichlet_value_[j] = voltage;
        }
        applyOuterBoundaryToSolution();
    }

    // A convenience: set *all* masked outer nodes to a constant (e.g. ground all).
    void setAllOuterBoundary(double voltage) {
        for (int j = 0; j < Nt_; ++j) {
            if (outer_dirichlet_mask_[j]) outer_dirichlet_value_[j] = voltage;
        }
        applyOuterBoundaryToSolution();
    }

    // Iterative solve for Laplace equation. Returns iteration count taken.
    // Stops when L_inf residual < tol or max_iters is reached.
    int run(int max_iters = 10'000, double tol = 1e-9) {
        // Keep boundary nodes fixed each iteration.
        // Center (r=0) uses symmetry: phi(0,theta) = average of first-ring neighbors.

        for (int iter = 0; iter < max_iters; ++iter) {
            double max_change = 0.0;

            // Enforce outer Dirichlet each sweep
            enforceOuterBC();

            // r = 0 (index i = 0): set to mean of i=1 ring to enforce regularity
            {
                double mean = 0.0;
                for (int j = 0; j < Nt_; ++j) mean += at(1, j);
                mean /= Nt_;
                for (int j = 0; j < Nt_; ++j) {
                    max_change = std::max(max_change, std::abs(at(0, j) - mean));
                    set(0, j, mean);
                }
            }

            // interior rings: 1 .. Nr-2 (since Nr-1 is boundary)
            for (int i = 1; i <= Nr_ - 2; ++i) {
                double r = i * dr_;
                double rpp = r + 0.5 * dr_;
                double rmm = r - 0.5 * dr_;
                double c_r = 1.0 / (dr_ * dr_);         // central radial coefficient (approximate)
                double c_th = 1.0 / (r * r * dth_ * dth_);

                for (int j = 0; j < Nt_; ++j) {
                    // neighbors
                    double phi_im = at(i - 1, j);
                    double phi_ip = at(i + 1, j);
                    double phi_jm = at(i, jm(j));
                    double phi_jp = at(i, jp(j));

                    // Conservative-ish FV/FD mix:
                    // (1/r) d/dr (r dphi/dr) ~ (rpp(phi_ip-phi_i) - rmm(phi_i-phi_im)) / (r * dr^2)
                    // => radial part: ((rpp * phi_ip + rmm * phi_im) - (rpp + rmm) * phi_i) / (r * dr^2)
                    // angular part: (phi_{j+1} - 2 phi_i + phi_{j-1}) / (r^2 dth^2)
                    double A = (rpp + rmm) / (r * dr_ * dr_) + 2.0 * c_th; // diagonal coeff
                    double rhs = (rpp * phi_ip + rmm * phi_im) / (r * dr_ * dr_)
                               + c_th * (phi_jp + phi_jm);

                    double phi_new = rhs / A;

                    if (useSOR_) {
                        double phi_old = at(i, j);
                        phi_new = phi_old + omega_ * (phi_new - phi_old);
                    }

                    max_change = std::max(max_change, std::abs(phi_new - at(i, j)));
                    set(i, j, phi_new);
                }
            }

            // re-enforce outer ring after interior update
            enforceOuterBC();

            if (max_change < tol) return iter + 1;
        }
        return max_iters;
    }

    // Query potential at (x,y) via bilinear interpolation in (r,theta) space.
    double getPotential(double x, double y) const {
        double r = std::hypot(x, y);
        if (r > R_) {
            // Outside: return nearest boundary value for robustness
            double th = std::atan2(y, x);
            int j = thetaToIndex(th);
            return atConst(Nr_ - 1, j);
        }
        // map to fractional indices
        double fr = r / dr_;
        int    i0 = std::clamp((int)std::floor(fr), 0, Nr_ - 2);
        double tr = fr - i0;

        double th = thetaWrap(std::atan2(y, x));
        double ft = th / dth_;
        int    j0 = (int)std::floor(ft) % Nt_;
        if (j0 < 0) j0 += Nt_;
        double tt = ft - std::floor(ft);

        // four neighbors (with periodic theta)
        int j1 = (j0 + 1) % Nt_;
        double p00 = atConst(i0,     j0);
        double p10 = atConst(i0 + 1, j0);
        double p01 = atConst(i0,     j1);
        double p11 = atConst(i0 + 1, j1);

        double p0 = p00 * (1 - tr) + p10 * tr;
        double p1 = p01 * (1 - tr) + p11 * tr;
        return p0 * (1 - tt) + p1 * tt;
    }

    // Utility: get raw grid (for debugging/printing)
    const std::vector<double>& rawGrid() const { return phi_; }
    int Nr() const { return Nr_; }
    int Ntheta() const { return Nt_; }
    double radius() const { return R_; }
    double dr() const { return dr_; }
    double dtheta() const { return dth_; }

private:
    // Grid & parameters
    double R_;
    int Nr_, Nt_;
    double dr_, dth_;
    bool useSOR_;
    double omega_;

    // Solution
    std::vector<double> phi_, phi_new_;

    // Outer boundary masks/values (only store outer ring, j=0..Nt-1)
    std::vector<bool>   outer_dirichlet_mask_;
    std::vector<double> outer_dirichlet_value_;

    // Bookkeeping per-electrode (stores theta indices on the boundary)
    std::vector<std::vector<int>> electrodes_;

    // Indexing helpers
    inline int idx(int i, int j) const { return i * Nt_ + (j % Nt_ + Nt_) % Nt_; }
    inline int jp(int j) const { return (j + 1) % Nt_; }
    inline int jm(int j) const { return (j - 1 + Nt_) % Nt_; }

    inline double& at(int i, int j) { return phi_[idx(i, j)]; }
    inline const double& atConst(int i, int j) const { return phi_[idx(i, j)]; }
    inline void set(int i, int j, double v) { phi_[idx(i, j)] = v; }

    inline int thetaToIndex(double theta) const {
        double t = thetaWrap(theta);
        int j = (int)std::round(t / dth_);
        j %= Nt_;
        if (j < 0) j += Nt_;
        return j;
    }
    inline double thetaWrap(double theta) const {
        double t = std::fmod(theta, 2.0*M_PI);
        return (t < 0.0) ? t + 2.0*M_PI : t;
        }

    /* // After changing masks/values, push them to the solution’s boundary
    void applyOuterBoundaryToSolution() {
        for (int j = 0; j < Nt_; ++j) {
            if (outer_dirichlet_mask_[j]) {
                set(Nr_ - 1, j, outer_dirichlet_value_[j]);
            }
        }
    }

    void enforceOuterBC() {
        for (int j = 0; j < Nt_; ++j) {
            if (outer_dirichlet_mask_[j]) {
                set(Nr_ - 1, j, outer_dirichlet_value_[j]);
            }
        }
    } */

    // After changing masks/values, push them to the solution’s boundary.
    // Dirichlet on masked indices; Neumann (zero-flux) elsewhere.
    void applyOuterBoundaryToSolution() {
        for (int j = 0; j < Nt_; ++j) {
            if (outer_dirichlet_mask_[j]) {
                set(Nr_ - 1, j, outer_dirichlet_value_[j]);             // Dirichlet
            } else {
                set(Nr_ - 1, j, at(Nr_ - 2, j));                        // Neumann: ∂φ/∂r=0
            }
        }
    }

    // Enforce boundary conditions every iteration sweep.
    void enforceOuterBC() {
        for (int j = 0; j < Nt_; ++j) {
            if (outer_dirichlet_mask_[j]) {
                set(Nr_ - 1, j, outer_dirichlet_value_[j]);             // Dirichlet
            } else {
                set(Nr_ - 1, j, at(Nr_ - 2, j));                        // Neumann: ∂φ/∂r=0
            }
        }
    }
    // Snapshot current mask indices for an "electrode id"
    std::vector<int> snapshotElectrodeIndices() const {
        std::vector<int> js;
        js.reserve(Nt_);
        for (int j = 0; j < Nt_; ++j)
            if (outer_dirichlet_mask_[j]) js.push_back(j);
        return js;
    }
};