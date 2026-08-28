// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko and David Llamas

#pragma once

#include <optional>
#include <string>

struct SimulationParameters {
    explicit SimulationParameters(double requiredHoppingDistance)
        : hoppingDistance(requiredHoppingDistance) {}

    double hoppingDistance;
    unsigned int totalNumSpins=10000;
    std::optional<unsigned int> gridSize;
    double hours=0.0;
    unsigned long maxSteps=1000000000000UL;
    std::optional<unsigned long> temperatureRoundSteps; // initial number of proposals between cooling updates
    double initialTemperature=20.0;
    double finalTemperature=0.1;
    int annealingSteps=1000;                             // number of cooling steps from initial to final temperature
    std::string initialConfiguration="random";
    int deltaOption=0;                                   // selects one of the two delta-function discretizations
    std::optional<unsigned long> randomSeed;
    bool overwriteExistingOutputs=false;
};

SimulationParameters parseSimulationParameters(int argc, char* argv[]);
