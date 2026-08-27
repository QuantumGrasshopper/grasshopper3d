#pragma once

#include <optional>
#include <string>

struct SimulationParameters {
    explicit SimulationParameters(double requiredHoppingDistance)
        : hoppingDistance(requiredHoppingDistance) {}

    double hoppingDistance;
    unsigned int totalNumSpins=5000;
    std::optional<unsigned int> gridSize;
    double hours=0.0;
    unsigned long maxSteps=1000000000000UL;
    std::optional<unsigned long> temperatureRoundSteps;
    double initialTemperature=25.0;
    double finalTemperature=0.1;
    int annealingSteps=1000;
    std::string initialConfiguration="random";
    int deltaOption=0;
    std::optional<unsigned long> randomSeed;
};

SimulationParameters parseSimulationParameters(int argc, char* argv[]);
