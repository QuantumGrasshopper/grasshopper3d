// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko

#include "doctest/doctest.h"
#include "setup.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace {

class ScopedTemporaryDirectory {
public:
    ScopedTemporaryDirectory()
        : originalDirectory_(std::filesystem::current_path())
        {
        const auto identifier=std::chrono::steady_clock::now().time_since_epoch().count();
        const std::string baseName="grasshopper3d-setup-tests-"+std::to_string(identifier);

        for(unsigned int attempt=0;attempt<100;attempt++)
            {
            const auto candidate=std::filesystem::temp_directory_path()
                               /(baseName+"-"+std::to_string(attempt));
            if(std::filesystem::create_directory(candidate))
                {
                temporaryDirectory_=candidate;
                std::filesystem::current_path(temporaryDirectory_);
                return;
                }
            }

        throw std::runtime_error("Could not create a temporary test directory");
        }

    ~ScopedTemporaryDirectory()
        {
        std::error_code error;
        std::filesystem::current_path(originalDirectory_,error);
        std::filesystem::remove_all(temporaryDirectory_,error);
        }

    ScopedTemporaryDirectory(const ScopedTemporaryDirectory&)=delete;
    ScopedTemporaryDirectory& operator=(const ScopedTemporaryDirectory&)=delete;

private:
    std::filesystem::path originalDirectory_;
    std::filesystem::path temporaryDirectory_;
};

void configureLoadTest()
    {
    totalNumSpins=3;
    cellSize=1.0;
    gridSize=2;
    gridVolume=8;
    tempScaling=1.0;
    deltaOption=0;
    }

void writeInitialConfiguration(const std::string& contents)
    {
    std::ofstream configuration("initconf.dat");
    REQUIRE(configuration.is_open());
    configuration << contents;
    }

void checkLoadRejected(const std::string& contents, const std::string& expectedText)
    {
    writeInitialConfiguration(contents);
    std::vector<unsigned char> grid(gridVolume);
    std::vector<int> spins(totalNumSpins);

    try
        {
        initLoad(grid.data(),spins.data());
        FAIL_CHECK("invalid initial configuration was accepted");
        }
    catch(const std::runtime_error& error)
        {
        CHECK(std::string(error.what()).find(expectedText)!=std::string::npos);
        }
    }

} // namespace

TEST_CASE("valid loaded configuration has exact occupancy")
    {
    configureLoadTest();
    ScopedTemporaryDirectory temporaryDirectory;
    writeInitialConfiguration("0\n3\n7\n");

    std::vector<unsigned char> grid(gridVolume,1);
    std::vector<int> spins(totalNumSpins,-1);
    initLoad(grid.data(),spins.data());

    CHECK(spins==std::vector<int>{0,3,7});
    CHECK(std::count(grid.begin(),grid.end(),static_cast<unsigned char>(1))
          ==static_cast<std::ptrdiff_t>(totalNumSpins));
    CHECK(grid[0]==1);
    CHECK(grid[3]==1);
    CHECK(grid[7]==1);
    }

TEST_CASE("configuration loading rejects invalid file contents")
    {
    configureLoadTest();
    ScopedTemporaryDirectory temporaryDirectory;

    SUBCASE("duplicate")
        {
        checkLoadRejected("0\n0\n7\n","Duplicate coordinate");
        }
    SUBCASE("out of range")
        {
        checkLoadRejected("0\n3\n8\n","outside the current grid");
        }
    SUBCASE("insufficient")
        {
        checkLoadRejected("0\n3\n","Invalid or insufficient data");
        }
    SUBCASE("extra")
        {
        checkLoadRejected("0\n3\n7\n6\n","more data than expected");
        }
    SUBCASE("malformed")
        {
        checkLoadRejected("0\nword\n7\n","Invalid or insufficient data");
        }
    }

TEST_CASE("configuration loading rejects a missing file")
    {
    configureLoadTest();
    ScopedTemporaryDirectory temporaryDirectory;
    std::vector<unsigned char> grid(gridVolume);
    std::vector<int> spins(totalNumSpins);

    CHECK_THROWS_WITH_AS(initLoad(grid.data(),spins.data()),
                         "Error: Cannot open initconf.dat for reading.",
                         std::runtime_error);
    }
