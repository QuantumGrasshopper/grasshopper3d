#include "interactions.hpp"
#include "utilities.hpp"

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;

bool isAround(double have, double comparewith)
	{
	/// Tests whether two physical lengths differ by at most two cell widths,
	/// which is the compact support of the smeared radial delta function.
	if(abs(have-comparewith)/cellSize<=2) return true;
	else return false;
	}

/// Grasshopper interaction contribution using the globally selected deltaOption
/// The inputs are the distance between two points and the interaction distance (in any order)
double contributionEnergy(double have, double comparewith)
	{
	double contribution=0;

	if(isAround(have,comparewith))
        {
        if(deltaOption==0) contribution=(1. + cos(PI*(have-comparewith)/cellSize/2.))/4.;
        else if(deltaOption==1)
            {
		double absdist=abs(have-comparewith)/cellSize;
            if(absdist<1) contribution=17./48.+sqrt(3.)*PI/108.+absdist/4.-absdist*absdist/4.+(1-2*absdist)*sqrt(1.+12*absdist*(1-absdist))/16.-sqrt(3.)*asin(sqrt(3.)*(2*absdist-1)/2.)/12.;
            else if( (absdist>=1)&&(absdist<2) ) contribution=55./48.-sqrt(3.)*PI/108.-13.*absdist/12.+absdist*absdist/4.+(2*absdist-3)*sqrt(36*absdist-23.-12*absdist*absdist)/48.+sqrt(3.)*asin(sqrt(3.)*(2*absdist-3)/2.)/36.;
            }
        else throw logic_error("Invalid delta function discretization option");
        }

	return contribution;
	}

int xcoord(int gridPoint)
    {
    // Flattened coordinates use index = z*gridSize^2 + y*gridSize + x.
    const int signedGridSize=static_cast<int>(gridSize);
    return gridPoint%signedGridSize;
    }

int ycoord(int gridPoint)
    {
    const int signedGridSize=static_cast<int>(gridSize);
    return gridPoint/signedGridSize-signedGridSize*zcoord(gridPoint);
    }

int zcoord(int gridPoint)
    {
    const int signedGridSize=static_cast<int>(gridSize);
    return gridPoint/signedGridSize/signedGridSize;
    }

int getGridPoint(int x, int y, int z)
    {
    const int signedGridSize=static_cast<int>(gridSize);
    return z*signedGridSize*signedGridSize + y*signedGridSize + x;
    }

tuple<double,double,double> findPosition(int gridPoint)
	{
    const int signedGridSize=static_cast<int>(gridSize);
    int z = zcoord(gridPoint);
    int y = gridPoint/signedGridSize-signedGridSize*z;
    int x = gridPoint - signedGridSize*signedGridSize*z - signedGridSize*y;

	tuple<double,double,double> thisPair(x*cellSize,y*cellSize,z*cellSize);
	return thisPair;
	}

double euclideanDistance(tuple<double,double,double> point1, tuple<double,double,double> point2)
	{
    double distance = pow(get<0>(point1)-get<0>(point2),2) + pow(get<1>(point1)-get<1>(point2),2) + pow(get<2>(point1)-get<2>(point2),2);
	return sqrt(distance);
	}

double euclideanDistance(tuple<int,int,int> point1, tuple<int,int,int> point2)
	{
    double distance = pow(get<0>(point1)-get<0>(point2),2) + pow(get<1>(point1)-get<1>(point2),2) + pow(get<2>(point1)-get<2>(point2),2);
	return cellSize*sqrt(distance);
	}

void validateInteractionTemplateReach(double distance)
    {
    // The centered template must contain every displacement with |r-distance| < 2*cellSize
    const unsigned int templateReach=(gridSize-1)/2;
    const double requiredReach=ceil(distance/cellSize)+1.;
    if(double(templateReach)<requiredReach)
        throw invalid_argument("Grid size is too small for the requested interaction-distance support.");
    }

GrasshopperInteractionTemplate buildInteractionTemplate(double distance)
    {
    validateInteractionTemplateReach(distance);

    // Build one template around the grid center. It stores only relative 3D
    // offsets with nonzero contributions, rather than a full translated table.
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
    // interactionGrid[i] = sum_j s_j w(i,j) for every grid site i,
    // including empty sites.
    vector<double> interactionGrid(gridVolume);
    const int signedGridSize=static_cast<int>(gridSize);

    // Translate the relative template to each site. Off-grid translations are
    // omitted, implementing open rather than periodic boundaries.
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

    // order is retained for floating-point reproducibility
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
    // Summing the cached field over occupied sites counts every occupied pair twice.
    return interaction/2.;
    }
