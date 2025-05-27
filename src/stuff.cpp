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