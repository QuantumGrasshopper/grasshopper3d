#include "interactions.hpp"
#include "utilities.hpp"

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;

void validateInteractionTemplateReach(double distance)
    {
    // The centered template must contain every displacement with
    // |r-distance| < 2*cellSize. Its shorter reach controls both even and odd grids.
    const unsigned int templateReach=(gridSize-1)/2;
    const double requiredReach=ceil(distance/cellSize)+1.;
    if(double(templateReach)<requiredReach)
        throw invalid_argument("Grid size is too small for the requested interaction-distance support.");
    }

GrasshopperInteractionTemplate buildInteractionTemplate(double distance)
    {
    validateInteractionTemplateReach(distance);

    GrasshopperInteractionTemplate interactionTemplate;
    const int signedGridSize=static_cast<int>(gridSize);
    const int centerCoordinate=signedGridSize/2;
    const tuple<int,int,int> centerPosition=
        make_tuple(centerCoordinate,centerCoordinate,centerCoordinate);

    for(int x=0;x<signedGridSize;x++)
        {
        for(int y=0;y<signedGridSize;y++)
            {
            for(int z=0;z<signedGridSize;z++)
                {
                const double contribution=contributionEnergy(
                    distance,euclideanDistance(centerPosition,make_tuple(x,y,z)));
                if(contribution > EPS)
                    {
                    RelativeNeighbor relativeNeighbor;
                    relativeNeighbor.dx=x-centerCoordinate;
                    relativeNeighbor.dy=y-centerCoordinate;
                    relativeNeighbor.dz=z-centerCoordinate;
                    relativeNeighbor.contribution=contribution;
                    interactionTemplate.push_back(relativeNeighbor);
                    }
                }
            }
        }

    return interactionTemplate;
    }

vector<double> buildGrasshopperInteractionGrid(
    const unsigned char grid[], const GrasshopperInteractionTemplate& interactionTemplate)
    {
    vector<double> interactionGrid(gridVolume);
    const int signedGridSize=static_cast<int>(gridSize);

    for(int siteX=0;siteX<signedGridSize;siteX++)
        {
        for(int siteY=0;siteY<signedGridSize;siteY++)
            {
            for(int siteZ=0;siteZ<signedGridSize;siteZ++)
                {
                const int site=getGridPoint(siteX,siteY,siteZ);
                const size_t siteIndex=static_cast<size_t>(site);
                interactionGrid[siteIndex]=0;
                for(size_t neighbor=0;neighbor<interactionTemplate.size();neighbor++)
                    {
                    const int neighborX=siteX+interactionTemplate[neighbor].dx;
                    const int neighborY=siteY+interactionTemplate[neighbor].dy;
                    const int neighborZ=siteZ+interactionTemplate[neighbor].dz;
                    if(neighborX >= 0 && neighborY >= 0 && neighborZ >= 0
                       && neighborX < signedGridSize && neighborY < signedGridSize
                       && neighborZ < signedGridSize)
                        {
                        if(grid[getGridPoint(neighborX,neighborY,neighborZ)]==true)
                            interactionGrid[siteIndex]+=interactionTemplate[neighbor].contribution;
                        }
                    }
                }
            }
        }

    return interactionGrid;
    }

double totalGrasshopperInteraction(
    const unsigned char grid[], const vector<double>& interactionGrid)
    {
    double interaction=0;
    const int signedGridSize=static_cast<int>(gridSize);

    // Preserve the original x/y/z summation order.
    for(int x=0;x<signedGridSize;x++)
        {
        for(int y=0;y<signedGridSize;y++)
            {
            for(int z=0;z<signedGridSize;z++)
                {
                const int site=getGridPoint(x,y,z);
                if(grid[site]==true)
                    interaction+=interactionGrid[static_cast<size_t>(site)];
                }
            }
        }

    return interaction/2.;
    }
