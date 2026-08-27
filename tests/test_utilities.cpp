#include "doctest/doctest.h"
#include "utilities.hpp"

#include <cmath>
#include <initializer_list>
#include <tuple>

namespace {

double symmetricKernelNormalization(int option)
    {
    cellSize=1.0;
    deltaOption=option;
    constexpr int intervals=200000;
    constexpr double upperLimit=2.0;
    const double step=upperLimit/intervals;
    double weightedSum=contributionEnergy(0.0,0.0)
                      +contributionEnergy(0.0,upperLimit);

    for(int interval=1;interval<intervals;interval++)
        {
        const double value=contributionEnergy(0.0,interval*step);
        weightedSum+=(interval%2==0 ? 2.0 : 4.0)*value;
        }

    return 2.0*step*weightedSum/3.0;
    }

} // namespace

TEST_CASE("3D coordinate conversion and distance geometry") {
    cellSize = 0.25;
    gridSize = 4;
    gridVolume = 64;

    const int gridPoint = getGridPoint(3, 2, 1);
    CHECK(gridPoint == 27);
    CHECK(xcoord(gridPoint) == 3);
    CHECK(ycoord(gridPoint) == 2);
    CHECK(zcoord(gridPoint) == 1);

    const auto position = findPosition(gridPoint);
    CHECK(std::get<0>(position) == doctest::Approx(0.75));
    CHECK(std::get<1>(position) == doctest::Approx(0.5));
    CHECK(std::get<2>(position) == doctest::Approx(0.25));

    CHECK(euclideanDistance(std::tuple<double,double,double>{0.0, 0.0, 0.0},
                            std::tuple<double,double,double>{3.0, 4.0, 12.0})
          == doctest::Approx(13.0));
    CHECK(euclideanDistance(std::tuple<int,int,int>{0, 0, 0},
                            std::tuple<int,int,int>{3, 4, 12})
          == doctest::Approx(3.25));
}

TEST_CASE("delta option zero has the expected compact support and weights") {
    cellSize = 0.25;
    deltaOption = 0;

    CHECK(isAround(0.0, 0.5));
    CHECK_FALSE(isAround(0.0, 0.500001));
    CHECK(contributionEnergy(0.0, 0.0) == doctest::Approx(0.5));
    CHECK(contributionEnergy(0.0, 0.25) == doctest::Approx(0.25));
    CHECK(std::abs(contributionEnergy(0.0, 0.5)) < 1e-15);
    CHECK(contributionEnergy(0.0, 0.500001) == 0.0);
    CHECK(contributionEnergy(0.25, 0.0) == contributionEnergy(0.0, 0.25));
}

TEST_CASE("delta option one evaluates both analytic branches and compact support") {
    cellSize = 1.0;
    deltaOption = 1;

    CHECK(contributionEnergy(0.0, 0.5)
          == doctest::Approx(0.46704998233983941).epsilon(1e-12));
    CHECK(contributionEnergy(0.0, 1.5)
          == doctest::Approx(0.03295001766016048).epsilon(1e-12));
    CHECK(contributionEnergy(0.0, 2.0) == 0.0);
    CHECK(contributionEnergy(0.0, 2.000001) == 0.0);
    CHECK(contributionEnergy(1.5, 0.0) == contributionEnergy(0.0, 1.5));
}

TEST_CASE("both radial kernels have unit symmetric normalization") {
    for(int option : {0,1})
        {
        CAPTURE(option);
        CHECK(symmetricKernelNormalization(option)
              ==doctest::Approx(1.0).epsilon(1e-10));
        }
}
