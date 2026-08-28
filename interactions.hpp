// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko and David Llamas

#pragma once

#include <tuple>
#include <vector>

bool isAround(double have, double comparewith);
double contributionEnergy(double have, double comparewith);
int xcoord(int gridPoint);
int ycoord(int gridPoint);
int zcoord(int gridPoint);
int getGridPoint(int x, int y, int z);
std::tuple<double,double,double> findPosition(int gridPoint);
double euclideanDistance(std::tuple<double,double,double> point1, std::tuple<double,double,double> point2);
double euclideanDistance(std::tuple<int,int,int> point1, std::tuple<int,int,int> point2);

struct RelativeNeighbor {
    int dx;
    int dy;
    int dz;
    double contribution;
};

using GrasshopperInteractionTemplate = std::vector<RelativeNeighbor>;

void validateInteractionTemplateReach(double distance);
GrasshopperInteractionTemplate buildInteractionTemplate(double distance);
std::vector<double> buildGrasshopperInteractionGrid(
    const unsigned char grid[], const GrasshopperInteractionTemplate& interactionTemplate);
double totalGrasshopperInteraction(
    const unsigned char grid[], const std::vector<double>& interactionGrid);
