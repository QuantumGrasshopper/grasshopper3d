#pragma once

#include <gsl/gsl_rng.h>
#include <tuple>

#define PI 3.14159265358979323846264338328 
#define EPS 1e-8

extern unsigned int totalNumSpins;
extern double cellSize;
extern unsigned int gridSize;
extern unsigned int gridVolume;
extern double tempScaling;
extern int deltaOption;

bool isAround(double have, double comparewith);
double contributionEnergy(double have, double comparewith);
int xcoord(int gridPoint);
int ycoord(int gridPoint);
int zcoord(int gridPoint);
int getGridPoint(int x, int y, int z);
std::tuple<double,double,double> findPosition(int gridPoint);
double euclideanDistance(std::tuple<double,double,double> point1, std::tuple<double,double,double> point2);
double euclideanDistance(std::tuple<int,int,int> point1, std::tuple<int,int,int> point2);
