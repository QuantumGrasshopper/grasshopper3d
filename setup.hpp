// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko and David Llamas

#pragma once

#include <gsl/gsl_rng.h>

#include <string>

void initLoad(unsigned char grid[], int spinArray[]);
void initRandom(unsigned char grid[], int spinArray[], gsl_rng *RNG);
void saveConfig(const int spinArray[], const std::string& filename);
void initialize(unsigned char grid[], int spinArray[], gsl_rng *RNG, std::string initconf);
