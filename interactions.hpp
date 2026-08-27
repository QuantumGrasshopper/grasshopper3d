#pragma once

#include <vector>

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
