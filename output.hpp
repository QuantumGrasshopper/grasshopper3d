// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko

#pragma once

#include <cstddef>
#include <fstream>
#include <limits>
#include <ostream>
#include <string>

constexpr int outputPrecision=std::numeric_limits<double>::max_digits10;
// Hard safety ceiling for optional trajectory snapshots; final configuration files are separate.
constexpr std::size_t maximumConfigurationFileBytes=100ULL*1024ULL*1024ULL;

void prepareOutputFiles(bool overwrite, bool preserveInitialConfiguration);
void setFullOutputPrecision(std::ostream& stream);
void openOutputFile(std::ofstream& stream, const std::string& filename);
void checkOutputStream(const std::ostream& stream,
                       const std::string& filename,
                       const char* operation);
void finishOutputFile(std::ofstream& stream, const std::string& filename);

class BoundedOutputFile {
public:
    BoundedOutputFile(const std::string& filename, std::size_t maximumBytes);
    bool writeIfFits(const std::string& data);
    bool limitReached() const;
    std::size_t bytesWritten() const;
    void finish();

private:
    std::ofstream file;
    std::string filename;
    std::size_t maximumBytes;
    std::size_t writtenBytes=0;
    bool reachedLimit=false;
};
