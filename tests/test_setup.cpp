// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko

#include "doctest/doctest.h"
#include "setup.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
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

void configureBallTest(unsigned int spins, unsigned int size)
    {
    totalNumSpins=spins;
    cellSize=1.0;
    gridSize=size;
    gridVolume=gridSize*gridSize*gridSize;
    tempScaling=1.0;
    deltaOption=0;
    }

uint64_t squaredDistanceFromCenter(unsigned int gridPoint)
    {
    const uint64_t planeSize=static_cast<uint64_t>(gridSize)*gridSize;
    const int64_t x=gridPoint%gridSize;
    const int64_t y=(gridPoint/gridSize)%gridSize;
    const int64_t z=gridPoint/planeSize;
    const int64_t doubledCenter=static_cast<int64_t>(gridSize)-1;
    const int64_t dx=2*x-doubledCenter;
    const int64_t dy=2*y-doubledCenter;
    const int64_t dz=2*z-doubledCenter;
    return static_cast<uint64_t>(dx*dx+dy*dy+dz*dz);
    }

void checkBallConfiguration(const std::vector<unsigned char>& grid,
                            const std::vector<int>& spins)
    {
    REQUIRE(grid.size()==gridVolume);
    REQUIRE(spins.size()==totalNumSpins);
    CHECK(std::is_sorted(spins.begin(),spins.end()));

    const std::set<int> uniqueSpins(spins.begin(),spins.end());
    CHECK(uniqueSpins.size()==totalNumSpins);
    CHECK(std::count(grid.begin(),grid.end(),static_cast<unsigned char>(1))
          ==static_cast<std::ptrdiff_t>(totalNumSpins));

    for(int spin : spins)
        {
        REQUIRE(spin>=0);
        REQUIRE(spin<static_cast<int>(gridVolume));
        CHECK(grid[static_cast<size_t>(spin)]==1);
        }
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

TEST_CASE("ball initialization selects the center and six axial neighbors on an odd grid")
    {
    configureBallTest(7,5);
    std::vector<unsigned char> grid(gridVolume);
    std::vector<int> spins(totalNumSpins);

    initBall(grid.data(),spins.data());

    checkBallConfiguration(grid,spins);
    CHECK(spins==std::vector<int>{37,57,61,62,63,67,87});
    }

TEST_CASE("ball initialization selects the central block on an even grid")
    {
    configureBallTest(8,4);
    std::vector<unsigned char> grid(gridVolume);
    std::vector<int> spins(totalNumSpins);

    initBall(grid.data(),spins.data());

    checkBallConfiguration(grid,spins);
    CHECK(spins==std::vector<int>{21,22,25,26,37,38,41,42});
    }

TEST_CASE("partial ball shells use a deterministic nearest-site tie break")
    {
    configureBallTest(4,5);
    std::vector<unsigned char> firstGrid(gridVolume);
    std::vector<unsigned char> secondGrid(gridVolume);
    std::vector<int> firstSpins(totalNumSpins);
    std::vector<int> secondSpins(totalNumSpins);

    initBall(firstGrid.data(),firstSpins.data());
    initBall(secondGrid.data(),secondSpins.data());

    checkBallConfiguration(firstGrid,firstSpins);
    checkBallConfiguration(secondGrid,secondSpins);
    CHECK(firstGrid==secondGrid);
    CHECK(firstSpins==secondSpins);
    CHECK(firstSpins==std::vector<int>{37,57,61,62});

    uint64_t maximumSelectedDistance=0;
    uint64_t minimumUnselectedDistance=std::numeric_limits<uint64_t>::max();
    for(unsigned int gridPoint=0;gridPoint<gridVolume;gridPoint++)
        {
        const uint64_t distance=squaredDistanceFromCenter(gridPoint);
        if(firstGrid[gridPoint])
            maximumSelectedDistance=std::max(maximumSelectedDistance,distance);
        else
            minimumUnselectedDistance=std::min(minimumUnselectedDistance,distance);
        }
    CHECK(maximumSelectedDistance<=minimumUnselectedDistance);
    }

TEST_CASE("the initialization path saves a generated ball without using the RNG")
    {
    configureBallTest(2,3);
    ScopedTemporaryDirectory temporaryDirectory;
    std::vector<unsigned char> grid(gridVolume);
    std::vector<int> spins(totalNumSpins);

    initialize(grid.data(),spins.data(),nullptr,"ball");

    checkBallConfiguration(grid,spins);
    std::ifstream savedConfiguration("initconf.dat");
    REQUIRE(savedConfiguration.is_open());
    std::vector<int> savedSpins;
    int savedSpin;
    while(savedConfiguration>>savedSpin) savedSpins.push_back(savedSpin);
    CHECK(savedSpins==spins);
    }
