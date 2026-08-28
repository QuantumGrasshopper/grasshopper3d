#include "setup.hpp"
#include "output.hpp"
#include "utilities.hpp"

#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>

using namespace std;

void initLoad(unsigned char grid[], int spinArray[])
{
	// Load the initial configuration file.
	ifstream initconfin("initconf.dat");
	if(!initconfin.is_open())
		throw runtime_error("Error: Cannot open initconf.dat for reading.");

	// Prepare the occupancy grid.
	for (unsigned int i = 0; i < gridVolume; i++)
	{
		grid[i] = false;
	}
	// Read exactly totalNumSpins unique flattened coordinates and validate
	// them against the current grid.
	for (unsigned int i = 0; i < totalNumSpins; i++)
	{
		long long coordinate;
		if(!(initconfin >> coordinate))
			throw runtime_error("Error: Invalid or insufficient data in initconf.dat.");
		if(coordinate<0 || coordinate>=static_cast<long long>(gridVolume))
			throw runtime_error("Error: Coordinate in initconf.dat is outside the current grid.");

		spinArray[i]=static_cast<int>(coordinate);
		if(grid[spinArray[i]]==true)
			throw runtime_error("Error: Duplicate coordinate in initconf.dat.");
		grid[spinArray[i]]=true;
	}

	// Reject trailing data after the expected coordinates.
	initconfin >> ws;
	if(initconfin.peek()!=char_traits<char>::eof())
		throw runtime_error("Error: initconf.dat contains more data than expected.");
}

void initRandom(unsigned char grid[], int spinArray[], gsl_rng *RNG)
{
	// Random initialization with exactly totalNumSpins distinct occupied cells.
	for (unsigned int i = 0; i < gridVolume; i++)
	{
		grid[i] = false;
	}
	unsigned long newSpinCoord;
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
		spinArray[spincounter] = static_cast<int>(newSpinCoord);
		spincounter++;
	}
}

void saveConfig(const int spinArray[], const string& filename)
{
	ofstream output;
	openOutputFile(output,filename);
	for (unsigned int i = 0; i < totalNumSpins; i++)
		output << spinArray[i] << '\n';
	finishOutputFile(output,filename);
}

void initialize(unsigned char grid[], int spinArray[], gsl_rng *RNG, string initconf)
{
	if (initconf == "random")
		initRandom(grid, spinArray, RNG);
	else if (initconf == "load")
		initLoad(grid, spinArray);
	else
		throw logic_error("Invalid initialization.");

	if (!(initconf == "load"))
	{
		saveConfig(spinArray,"initconf.dat");
	}
}
