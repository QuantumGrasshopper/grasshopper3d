// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko

#include "output.hpp"

#include <array>
#include <filesystem>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

using namespace std;

namespace {

const array<const char*,7> outputFiles{
    "result.dat",
    "energies.dat",
    "temperatures.dat",
    "config.dat",
    "initconf.dat",
    "finconf.dat",
    "bestconf.dat"
};

} // namespace

void prepareOutputFiles(bool overwrite, bool preserveInitialConfiguration)
    {
    vector<string> existingOutputFiles;

    // First inspect the complete output set without changing anything.
    for(const char* filename : outputFiles)
        {
        if(preserveInitialConfiguration && string(filename)=="initconf.dat") continue;

        error_code error;
        const filesystem::file_status status=filesystem::symlink_status(filename,error);
        if(error)
            {
            if(error==errc::no_such_file_or_directory) continue;
            throw runtime_error("Cannot inspect output artifact "+string(filename)
                                +": "+error.message());
            }
        if(!filesystem::exists(status)) continue;
        if(filesystem::is_directory(status))
            throw runtime_error("Output artifact is a directory: "+string(filename));
        if(!overwrite)
            throw runtime_error("Output artifact already exists: "+string(filename)
                                +". Use -overwrite 1 to replace existing outputs.");

        existingOutputFiles.push_back(filename);
        }

    // Only start removing files after the complete preflight succeeded.
    for(const string& filename : existingOutputFiles)
        {
        error_code error;
        const bool removed=filesystem::remove(filename,error);
        if(error || !removed)
            {
            const string detail=error ? error.message() : "file was not removed";
            throw runtime_error("Cannot remove output artifact "+filename+": "+detail);
            }
        }
    }

void setFullOutputPrecision(ostream& stream)
    {
    stream << setprecision(outputPrecision);
    }

void openOutputFile(ofstream& stream, const string& filename)
    {
    stream.open(filename);
    if(!stream.is_open())
        throw runtime_error("Failed to open output file "+filename+".");
    setFullOutputPrecision(stream);
    }

void checkOutputStream(const ostream& stream,
                       const string& filename,
                       const char* operation)
    {
    if(!stream)
        throw runtime_error("Failed to "+string(operation)
                            +" output file "+filename+".");
    }

void finishOutputFile(ofstream& stream, const string& filename)
    {
    checkOutputStream(stream,filename,"write");
    stream.flush();
    checkOutputStream(stream,filename,"flush");
    stream.close();
    checkOutputStream(stream,filename,"close");
    }

BoundedOutputFile::BoundedOutputFile(const string& filename,
                                     size_t maximumBytes)
    : filename(filename), maximumBytes(maximumBytes)
    {
    openOutputFile(file,filename);
    }

bool BoundedOutputFile::writeIfFits(const string& data)
    {
    if(reachedLimit) return false;
    if(data.size()>maximumBytes-writtenBytes)
        {
        reachedLimit=true;
        return false;
        }

    file << data;
    checkOutputStream(file,filename,"write");
    writtenBytes+=data.size();
    return true;
    }

bool BoundedOutputFile::limitReached() const
    {
    return reachedLimit;
    }

size_t BoundedOutputFile::bytesWritten() const
    {
    return writtenBytes;
    }

void BoundedOutputFile::finish()
    {
    finishOutputFile(file,filename);
    }
