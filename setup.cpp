#include "setup.hpp"

using namespace std;

void initLoad(bool grid[], int spinArray[])
{
	ifstream initconfin("initconf.dat");
	for (unsigned int i = 0; i < gridVolume; i++)
	{
		grid[i] = false;
	}
	for (unsigned int i = 0; i < totalNumSpins; i++)
	{
		initconfin >> spinArray[i];
		grid[spinArray[i]] = true;
	}
}

void initRandom(bool grid[], int spinArray[], gsl_rng *RNG)
{
	for (unsigned int i = 0; i < gridVolume; i++)
	{
		grid[i] = false;
	}
	int newSpinCoord;
	unsigned int spincounter = 0;
	while (spincounter < totalNumSpins)
	{
		bool create = true;
		while (create == true)
		{
			newSpinCoord = gsl_rng_uniform_int(RNG, gridVolume);
			create = grid[newSpinCoord];
		}
		grid[newSpinCoord] = true;
		spinArray[spincounter] = newSpinCoord;
		spincounter++;
	}
}

void saveConfig(int *spinArray, ofstream &filename)
{
	for (unsigned int i = 0; i < totalNumSpins; i++)
		filename << spinArray[i] << endl;
}

void initialize(bool grid[], int spinArray[], gsl_rng *RNG, string initconf)
{
	if (initconf == "random")
		initRandom(grid, spinArray, RNG);
	else if (initconf == "load")
		initLoad(grid, spinArray);
	else
		throw logic_error("Invalid initialization.");

	if (!(initconf == "load"))
	{
		ofstream initconfout("initconf.dat");
		saveConfig(spinArray, initconfout);
	}
}

void initialize_ball(bool grid[], int spinArray[])
{
	radius = cbrt((3 / 4) / PI);
	tuple<double, double, double> startposition = findPosition(((gridSize * gridSize) * gridSize / 2) + (gridSize * gridSize / 2) + (gridSize / 2));
	int countspins = 0;
	for (unsigned int i; i < gridVolume; i++)
	{
		grid[i] = false;
		if (countspins < totalNumSpins)
		{
			tuple<double, double, double> thispoint = findPosition(i);
			if (euclideanDistance(thispoint) < radius)
			{
				grid[i] = true;
				spinArray[countspins] = i;
				countspins++;
			}
		}
	}
	for (int i = 0; i < gridVolume; i++)
	{
		if ((countspins < totalNumSpins) && (grid[i] == false))
		{
			tuple<double, double, double> thisPoint = findPosition(i);
			if (euclideanDistance(startpoint, thisPoint) < radius + cellsize)
			{
				grid[i] = true;
				spinArray[countspins] = i;
				countspins++;
			}
		}
	}

	if (countspins != totalNumSpins)
	{
		result << "ERROR in countspins" << endl;
		break;
	}
}