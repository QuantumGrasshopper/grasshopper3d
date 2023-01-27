#include "common3d.hpp"

using namespace std;

bool isAround(double have, double comparewith)
	{	
	if(abs(have-comparewith)/cellsize<=2) return true;
	else return false;
	}
	
double contributionEnergy(double have, double comparewith)
	{
	double contribution=0;
	if(isAround(have,comparewith)) contribution=(1. + cos(PI*(have-comparewith)/cellsize/2.))/4.;
	return contribution;
	} 
	
int xcoord(int gridPoint)
    {
    return gridPoint%gridSize;
    }
    
int ycoord(int gridPoint)
    {
    //return ((gridPoint-xcoord(gridPoint))/gridSize)%gridSize;
    return gridPoint/gridSize-gridSize*zcoord(gridPoint);
    }
    
int zcoord(int gridPoint)
    {
    //return (gridPoint - xcoord(gridPoint) - gridSize*ycoord(gridPoint))/gridSize/gridSize;
    return gridPoint/gridSize/gridSize;
    }

tuple<double,double,double> findPosition(int gridPoint)
	{
    int z = zcoord(gridPoint);
    int y = gridPoint/gridSize-gridSize*z;
    int x = gridPoint - gridSize*gridSize*z - gridSize*y;
        
	tuple<double,double,double> thisPair(x*cellsize,y*cellsize,z*cellsize);
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
	return cellsize*sqrt(distance);
	}
	
