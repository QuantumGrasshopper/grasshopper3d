// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko

#include "utilities.hpp"
#include "setup.hpp"
#include "annealing.hpp"
#include "interactions.hpp"
#include "parameters.hpp"
#include "output.hpp"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

unsigned int totalNumSpins;
double cellSize;
unsigned int gridSize;
unsigned int gridVolume;
double tempScaling;
int deltaOption;

using namespace std; 

namespace {

unsigned int automaticGridSize(unsigned int numberSpins, double distance)
    {
	// If no explicit gridSize is supplied, choose a generous default box to reduce effects of artificial open boundaries
    const double bulkHalfSize=trunc(pow(double(numberSpins),1./3.)+EPS);
    const double interactionHalfSize=trunc(2*distance/cellSize+EPS);
    const double automaticSize=2*bulkHalfSize+2*interactionHalfSize;

    if(!isfinite(automaticSize) || automaticSize<1
       || automaticSize>numeric_limits<unsigned int>::max())
        throw invalid_argument("Automatically derived grid size is outside the supported range.");

    return static_cast<unsigned int>(automaticSize);
    }

unsigned int checkedGridVolume(unsigned int inputGridSize)
    {
    if(inputGridSize==0)
        throw invalid_argument("Grid size must be positive.");

    // Compute gridSize^3 in a wide type before narrowing, and require every
    // valid flattened coordinate to remain representable by the int-based helpers.
    const uint64_t side=inputGridSize;
    const uint64_t maximumVolume=uint64_t(numeric_limits<int>::max())+1;

    if(side>maximumVolume/side)
        throw invalid_argument("Grid size is outside the supported flattened-index range.");
    const uint64_t square=side*side;
    if(square>maximumVolume/side)
        throw invalid_argument("Grid size is outside the supported flattened-index range.");
    const uint64_t volume=square*side;

    if(volume>vector<double>().max_size())
        throw invalid_argument("Grid volume is outside the supported storage range.");

    return static_cast<unsigned int>(volume);
    }

}

int main(int inputN,char *inputV[]) {
    try {
    
    // SETUP -------------------------------------------------------------------------------------

	const SimulationParameters parameters=parseSimulationParameters(inputN,inputV);
	double d=parameters.hoppingDistance;					//grasshopper hopping distance
	double maxtime=60*60*parameters.hours*1000;
	long unsigned int maxsteps=parameters.maxSteps;
	totalNumSpins=parameters.totalNumSpins;
	long unsigned int temproundsteps=
		parameters.temperatureRoundSteps.value_or(totalNumSpins);
	double temperature=parameters.initialTemperature;
	double finaltemperature=parameters.finalTemperature;
	int numberannealingsteps=parameters.annealingSteps;
	string initconf=parameters.initialConfiguration;
    deltaOption=parameters.deltaOption;

	tempScaling=pow((finaltemperature/temperature),1./double(numberannealingsteps));
	int outputconfigbeforetherm=numberannealingsteps/100;
	if(outputconfigbeforetherm<1) outputconfigbeforetherm=1;
	// Snapshot scheduling for config.dat; BoundedOutputFile enforces a separate byte-size cap.
	int annealingcounter=0; int maxoutputconfigs=200;
	
	cellSize = pow(1/double(totalNumSpins),1./3.);
	if(parameters.gridSize.has_value()) gridSize=*parameters.gridSize;
	else gridSize=automaticGridSize(totalNumSpins,d);
	gridVolume = checkedGridVolume(gridSize);
	if(totalNumSpins==0 || totalNumSpins>=gridVolume)
		throw invalid_argument("Number of spins must satisfy 0 < N < grid volume.");
	validateInteractionTemplateReach(d);
	const int signedGridSize=static_cast<int>(gridSize);
    // Raw energy counts each occupied pair once. The corresponding
	// grasshopper success probability is P = E/(2*pi*d^2*N^(5/3)).
    double probabilityNormFactor = 1/2./PI/d/d/pow(double(totalNumSpins),5./3.);

	// Use the system clock if no random seed was supplied; the selected seed is always output.
	long unsigned int seed;
	if(parameters.randomSeed.has_value()) seed=*parameters.randomSeed;
    else seed=static_cast<unsigned long>(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count());

	vector<unsigned char> grid(gridVolume);			//true if spin=1 at this grid point
	vector<int> spinArray(totalNumSpins);			//grid point where any spin is
	// Validate a loaded configuration before existing output artifacts can be removed.
	if(initconf=="load")
		initLoad(grid.data(),spinArray.data());

	prepareOutputFiles(parameters.overwriteExistingOutputs,initconf=="load");

    ofstream result;
	openOutputFile(result,"result.dat");
    
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	unsigned long temproundcounter=0;
	double accratio;
	unsigned long counter=0;
	unsigned long accepted=0;
	unsigned long accepted_current=0;

    result << "3D Grasshopper with Simulated Annealing, Euclidean metric" << endl;
	result << endl;
	result << "Total number of spins: " << totalNumSpins << endl;
	result << "Hopping distance: " << d << endl;
	result << "Size of grid: " << gridSize << endl;
	result << "Size of cell: " << cellSize << endl;
	result << endl;
	result << "Random seed: " << seed << endl;
    result << "Option for delta-function discretization: " << deltaOption << endl;
	result << "Initial temperature: " << temperature << endl;
	result << "Final temperature: " << finaltemperature << endl;
	result << "Temperature scaling factor: " << tempScaling << endl;
	result << "Number of annealing steps: " << numberannealingsteps << endl;
	result << "Initial number of steps before temperature scaling: " << temproundsteps << endl;
	result << endl;
	checkOutputStream(result,"result.dat","write");

    auto begin = chrono::high_resolution_clock::now();
    
    // CONSTRUCT NEIGHBOR TEMPLATE ---------------------------------------------------------------------------  
    
	GrasshopperInteractionTemplate interactionTemplate=buildInteractionTemplate(d);
		
    auto now = chrono::high_resolution_clock::now();
    auto timeDiff = chrono::duration_cast<chrono::milliseconds>(now-begin).count();
	result << "Time to construct neighbors list: " << timeDiff << "ms" << endl;
	checkOutputStream(result,"result.dat","write");
    
    // INITIAL SPIN CONFIGURATION ------------------------------------------------------------------------

	vector<int> noSpinArray(gridVolume-totalNumSpins);	 //complementary to above: grid point where no spin is

	if(initconf!="load")
		initialize(grid.data(),spinArray.data(),RNG,initconf);
    
    size_t noSpinCounter=0;
    for(int i=0;i<signedGridSize;i++)
        {
        for(int j=0;j<signedGridSize;j++)
            {
            for(int n=0;n<signedGridSize;n++)
                {
                const int gridPoint=getGridPoint(i,j,n);
                if(grid[static_cast<size_t>(gridPoint)]==false)
                    {
                    if(noSpinCounter>=noSpinArray.size())
                        throw logic_error("Initialization produced an incorrect number of empty grid cells.");
                    noSpinArray[noSpinCounter]=gridPoint;
                    noSpinCounter++;
                    }
                }
            }
        }
	if(noSpinCounter!=noSpinArray.size())
		throw logic_error("Initialization produced an incorrect number of empty grid cells.");

	// energyGrid[i] is the interaction contribution of a hypothetical spin at i
	// with the current occupied set, whether site i is occupied or empty.
	vector<double> energyGrid=buildGrasshopperInteractionGrid(grid.data(),interactionTemplate);
	double energy=totalGrasshopperInteraction(grid.data(),energyGrid);

    ofstream energies;
	openOutputFile(energies,"energies.dat");
	ofstream temperatures;
	openOutputFile(temperatures,"temperatures.dat");
	energies << energy*probabilityNormFactor << '\n';
	checkOutputStream(energies,"energies.dat","write");

	BoundedOutputFile configuration("config.dat",maximumConfigurationFileBytes);
	// Each config.dat row contains N flattened coordinates followed by the grasshopper success probability.
	auto buildConfigurationSnapshot=[&]()
		{
		ostringstream buffer;
		setFullOutputPrecision(buffer);
		for(unsigned int i=0;i<totalNumSpins;i++) buffer << spinArray[i] << " ";
		buffer << energy*probabilityNormFactor << '\n';
		return buffer.str();
		};
	bool configurationOutputEnabled=
		configuration.writeIfFits(buildConfigurationSnapshot());
	
	vector<int> bestSpinArray(totalNumSpins);	//the overall best spin array during the whole run
	for(unsigned int i=0;i<totalNumSpins;i++) bestSpinArray[i]=spinArray[i];
	double bestenergy=energy;
		
    // MAIN LOOP ------------------------------------------------------------------------------------------
		
    while( (static_cast<double>(timeDiff)<maxtime) && (counter<maxsteps) )
        {
		counter++; temproundcounter++;
		
        // MC update
        // Choose one filled site and one empty site.
		size_t destroy=static_cast<size_t>(gsl_rng_uniform_int (RNG, totalNumSpins));
		int oldSpinCoord=spinArray[destroy];
		size_t create=static_cast<size_t>(gsl_rng_uniform_int (RNG, gridVolume-totalNumSpins));
		int newSpinCoord=noSpinArray[create];
            
        double energyDifference=energyGrid[static_cast<size_t>(newSpinCoord)]
            -energyGrid[static_cast<size_t>(oldSpinCoord)];
		// The cached field at the new site may include its interaction with the
		// still-occupied old site; subtract that pair because the old site will be emptied.
        if(isAround(d,euclideanDistance(findPosition(newSpinCoord),findPosition(oldSpinCoord))))//NOTE more efficient to keep this check explicit
            {
            energyDifference -= contributionEnergy(d,euclideanDistance( findPosition(newSpinCoord),findPosition(oldSpinCoord) ));
            }

		bool accept;
		if(energyDifference>=0) accept=true;
		else accept=acceptreject(energyDecreaseProbDistr(energyDifference,temperature),RNG);
		
		if(accept==true)
			{       
			// Update the occupied and empty site lists.
            grid[static_cast<size_t>(oldSpinCoord)]=false;
            grid[static_cast<size_t>(newSpinCoord)]=true;
			spinArray[destroy]=newSpinCoord; noSpinArray[create]=oldSpinCoord;
			energy+=energyDifference;
        
            // Update the cached field: remove the old site's contribution from
            // every affected value, then add the new site's contribution.
            for(size_t j=0;j<interactionTemplate.size();j++)
                {
                int x = interactionTemplate[j].dx+xcoord(oldSpinCoord);
                int y = interactionTemplate[j].dy+ycoord(oldSpinCoord);
                int z = interactionTemplate[j].dz+zcoord(oldSpinCoord);
                if(x >= 0 && y >= 0 && z >=0
                   && x < signedGridSize && y < signedGridSize && z < signedGridSize)
                    {
                    energyGrid[static_cast<size_t>(getGridPoint(x,y,z))] -= interactionTemplate[j].contribution;
                    }
                x = interactionTemplate[j].dx+xcoord(newSpinCoord);
                y = interactionTemplate[j].dy+ycoord(newSpinCoord);
                z = interactionTemplate[j].dz+zcoord(newSpinCoord);
                if(x >= 0 && y >= 0 && z >=0
                   && x < signedGridSize && y < signedGridSize && z < signedGridSize)
                    {
                    energyGrid[static_cast<size_t>(getGridPoint(x,y,z))] += interactionTemplate[j].contribution;
                    }
                }
                
			accepted_current++;
			// Keep track of the optimal configuration and its pair energy.
			if(energy>bestenergy)
				{
                bestenergy=energy; 
                for(unsigned int i=0;i<totalNumSpins;i++) bestSpinArray[i]=spinArray[i];
                }
			}
        
        // cooling step
		if(temproundcounter==temproundsteps)
			{
			if(temperature>finaltemperature) 
				{
				temperature=temperatureDecrease(temperature);
				if(configurationOutputEnabled
				   && annealingcounter%outputconfigbeforetherm==0)
					{
					configurationOutputEnabled=
						configuration.writeIfFits(buildConfigurationSnapshot());
					}
				annealingcounter++;
				}
			else if(annealingcounter<maxoutputconfigs)
				{
				annealingcounter++;
				if(configurationOutputEnabled)
					configurationOutputEnabled=
						configuration.writeIfFits(buildConfigurationSnapshot());
				}
			accratio=static_cast<double>(accepted_current)/double(temproundcounter);
			temperatures << counter << '\t' << temperature << '\t' << accratio << '\n';
			checkOutputStream(temperatures,"temperatures.dat","write");
			energies << energy*probabilityNormFactor << '\n';
			checkOutputStream(energies,"energies.dat","write");
			accepted+=accepted_current; accepted_current=0;
			temproundcounter=0;
			temproundsteps=stepIncrease(temproundsteps);
			}
		
		now = chrono::high_resolution_clock::now();
		timeDiff = chrono::duration_cast<chrono::milliseconds>(now-begin).count();
		}
		
    // WRAP UP --------------------------------------------------------------------------------------------

	const unsigned long totalAccepted=accepted+accepted_current;
	const double averageAcceptanceRatio=
		counter>0 ? static_cast<double>(totalAccepted)/double(counter) : 0.0;
    
    result << endl;
	result << "Simulation took " << static_cast<double>(timeDiff)/60./1000 << " minutes" << endl;
	result << "Finished after " << counter << " steps" << endl;
	result << "Final temperature: " << temperature << endl;
	result << "Average acceptance ratio: " << averageAcceptanceRatio << endl;
	result << endl;
	result << "final energy: " << energy << endl;
	result << "best energy: " << bestenergy << endl;
    result << "final probability: " << energy*probabilityNormFactor << endl;
	result << "best probability: " << bestenergy*probabilityNormFactor << endl;
	result << endl;
	checkOutputStream(result,"result.dat","write");

	finishOutputFile(energies,"energies.dat");
	finishOutputFile(temperatures,"temperatures.dat");
	configuration.finish();
	finishOutputFile(result,"result.dat");

    saveConfig(spinArray.data(),"finconf.dat");
    saveConfig(bestSpinArray.data(),"bestconf.dat");
    
    return 0;
    }
    catch (const exception& error) {
        cerr << "Error: " << error.what() << endl;
        return 1;
    }
}
