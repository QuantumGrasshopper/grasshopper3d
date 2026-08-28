#include "utilities.hpp"

#include <cmath>
#include <stdexcept>

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
	
