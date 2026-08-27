#include "utilities.hpp"
#include "setup.hpp"
#include "annealing.hpp"

#include <cstdint>
#include <limits>

unsigned int totalNumSpins;
double cellSize;
unsigned int gridSize;
unsigned int gridVolume;
double tempScaling;
int deltaOption;

using namespace std; 

struct gridcell {
    int x;
    int y;
    int z;
    double energy;    
};

namespace {

unsigned int automaticGridSize(unsigned int numberSpins, double distance)
    {
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

void validateInteractionTemplateReach(double distance)
    {
    // The centered template must contain every displacement with
    // |r-distance| < 2*cellSize. Its shorter reach controls both even and odd grids.
    const unsigned int templateReach=(gridSize-1)/2;
    const double requiredReach=ceil(distance/cellSize)+1.;
    if(double(templateReach)<requiredReach)
        throw invalid_argument("Grid size is too small for the requested interaction-distance support.");
    }

}

int main(int inputN,char *inputV[]) {
    
    // SETUP -------------------------------------------------------------------------------------
    
    ofstream result("result.dat");
	
	double d=get_option(inputN,inputV,"d");					//grasshopper hopping distance
	double maxtime=get_option(inputN,inputV,"hours");
	long unsigned int maxsteps=get_option(inputN,inputV,"steps");
	long unsigned int temproundsteps=get_option(inputN,inputV,"tempsteps");
	double temperature=get_option(inputN,inputV,"inittemp");
	double finaltemperature=get_option(inputN,inputV,"fintemp");
	int numberannealingsteps=get_option(inputN,inputV,"annealsteps");
	totalNumSpins=get_option(inputN,inputV,"N");
	gridSize=get_option(inputN,inputV,"gridsize");
	long unsigned int seed=get_option(inputN,inputV,"randomseed");
	string initconf=get_string_option(inputN,inputV,"initconf");
    deltaOption=get_option(inputN,inputV,"delta");
	
	maxtime=60*60*maxtime*1000;
	if(totalNumSpins<100) totalNumSpins=5000;
	if(maxsteps==0) maxsteps=1e12;
	if(temproundsteps>maxsteps) temproundsteps=int(maxsteps/1000.);
	if(temproundsteps<10) temproundsteps=totalNumSpins;
	if(temperature<EPS) temperature=25.;
	if(finaltemperature<EPS) finaltemperature=0.1;
	if(numberannealingsteps<EPS) numberannealingsteps=1000;
	tempScaling=pow((finaltemperature/temperature),1./double(numberannealingsteps));
	int outputconfigbeforetherm=numberannealingsteps/100; int annealingcounter=0; int maxoutputconfigs=200;	//for output config dat
	
	cellSize = pow(1/double(totalNumSpins),1./3.);
	if(gridSize<10) gridSize=automaticGridSize(totalNumSpins,d);
	gridVolume = checkedGridVolume(gridSize);
	if(totalNumSpins==0 || totalNumSpins>=gridVolume)
		throw invalid_argument("Number of spins must satisfy 0 < N < grid volume.");
	validateInteractionTemplateReach(d);
    // one factor of 1/2 is already taken care of by avoiding double counting
    double probabilityNormFactor = 1/2./PI/d/d/pow(double(totalNumSpins),5./3.);
    
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
    if(seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	gsl_rng_set (RNG, seed);
	
	unsigned int temproundcounter=0;
	double accratio;
    long unsigned int counter=0; long accepted=0; long accepted_current=0;

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

    auto begin = chrono::high_resolution_clock::now();
    
    // CONSTRUCT NEIGHBOR TEMPLATE ---------------------------------------------------------------------------  
    
    vector< gridcell > dNeighbourTemplate;
    
    double thisEnergyContribution;
    const int centerCoordinate=static_cast<int>(gridSize/2);
    tuple<int,int,int> centerPosition=make_tuple(centerCoordinate,centerCoordinate,centerCoordinate);
    for(unsigned int i=0;i<gridSize;i++) 
        {
        for(unsigned int j=0;j<gridSize;j++)
            {
            for(unsigned int n=0;n<gridSize;n++)
                {
                thisEnergyContribution=contributionEnergy(d,euclideanDistance(centerPosition,make_tuple(i,j,n)));
                if(thisEnergyContribution > EPS)
                    {
                    gridcell this_cell;
                    this_cell.x = static_cast<int>(i)-centerCoordinate;
                    this_cell.y = static_cast<int>(j)-centerCoordinate;
                    this_cell.z = static_cast<int>(n)-centerCoordinate;
                    this_cell.energy = thisEnergyContribution;
                    dNeighbourTemplate.push_back(this_cell);
                    }
                }
            }
        }
		
    auto now = chrono::high_resolution_clock::now();
    auto timeDiff = chrono::duration_cast<chrono::milliseconds>(now-begin).count();
	result << "Time to construct neighbors list: " << timeDiff << "ms" << endl;
    
    // INITIAL SPIN CONFIGURATION ------------------------------------------------------------------------
    
    vector<unsigned char> grid(gridVolume);			//true if spin=1 at this grid point
    vector<double> energyGrid(gridVolume);
	vector<int> spinArray(totalNumSpins);				 //grid point where any spin is
	vector<int> noSpinArray(gridVolume-totalNumSpins);	 //complementary to above: grid point where no spin is

    initialize(grid.data(), spinArray.data(), RNG, initconf);
    
    unsigned int noSpinCounter=0;
    double energy = 0;
    for(unsigned int i=0;i<gridSize;i++) 
        {
        for(unsigned int j=0;j<gridSize;j++)
            {
            for(unsigned int n=0;n<gridSize;n++)
                {
                const int gridPoint=getGridPoint(i,j,n);
                if(grid[gridPoint]==false) {noSpinArray[noSpinCounter]=gridPoint; noSpinCounter++;}
                // Fill the energy grid
                energyGrid[gridPoint]=0;
                for(unsigned int k=0;k<dNeighbourTemplate.size();k++)
                    {
                    int x = dNeighbourTemplate[k].x  +i;
                    int y = dNeighbourTemplate[k].y  +j;
                    int z = dNeighbourTemplate[k].z  +n;
                    if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                        {
                        if(grid[getGridPoint(x,y,z)]==true)
                            {energyGrid[gridPoint]+=dNeighbourTemplate[k].energy;}
                        }
                    }
                if(grid[gridPoint]==true) energy += energyGrid[gridPoint];
                }
            }
        }
    energy = energy/2.;
		
    ofstream energies("energies.dat");
	ofstream temperatures("temperatures.dat");
	energies << energy << endl;
	ofstream configuration("config.dat");
	for(unsigned int i=0;i<totalNumSpins;i++) configuration << spinArray[i] << " ";
	configuration << energy << endl;
	
	vector<int> bestSpinArray(totalNumSpins);	//the overall best spin array during the whole run
	for(unsigned int i=0;i<totalNumSpins;i++) bestSpinArray[i]=spinArray[i];
	double bestenergy=energy;
		
    // MAIN LOOP ------------------------------------------------------------------------------------------
		
    while( (timeDiff<maxtime) && (counter<maxsteps) )
        {
		counter++; temproundcounter++;
		
        // MC update
		int destroy=gsl_rng_uniform_int (RNG, totalNumSpins);
		int oldSpinCoord=spinArray[destroy];
		int create=gsl_rng_uniform_int (RNG, gridVolume-totalNumSpins);
		int newSpinCoord=noSpinArray[create];
            
        double energyDifference=energyGrid[newSpinCoord]-energyGrid[oldSpinCoord];
        if(isAround(d,euclideanDistance(findPosition(newSpinCoord),findPosition(oldSpinCoord))))//NOTE more efficient to keep this check explicit
            {
            energyDifference -= contributionEnergy(d,euclideanDistance( findPosition(newSpinCoord),findPosition(oldSpinCoord) ));
            }

		bool accept;
		if(energyDifference>=0) accept=true;
		else accept=acceptreject(energyDecreaseProbDistr(energyDifference,temperature),RNG);
		
		if(accept==true)
			{       
            grid[oldSpinCoord]=false; 
            grid[newSpinCoord]=true;        
			spinArray[destroy]=newSpinCoord; noSpinArray[create]=oldSpinCoord;
			energy+=energyDifference;
        
            //update energy grid                
            for(unsigned int j=0;j<dNeighbourTemplate.size();j++)
                {
                int x = dNeighbourTemplate[j].x  +xcoord(oldSpinCoord);
                int y = dNeighbourTemplate[j].y  +ycoord(oldSpinCoord);
                int z = dNeighbourTemplate[j].z  +zcoord(oldSpinCoord);
                if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                    {
                    energyGrid[getGridPoint(x,y,z)] -= dNeighbourTemplate[j].energy;
                    }
                x = dNeighbourTemplate[j].x  +xcoord(newSpinCoord);
                y = dNeighbourTemplate[j].y  +ycoord(newSpinCoord);
                z = dNeighbourTemplate[j].z  +zcoord(newSpinCoord);
                if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                    {
                    energyGrid[getGridPoint(x,y,z)] += dNeighbourTemplate[j].energy;
                    }
                }
                
			accepted_current++;
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
				if(annealingcounter%outputconfigbeforetherm==0) 
					{
                    for(unsigned int i=0;i<totalNumSpins;i++) configuration << spinArray[i] << " "; 
                    configuration << energy << endl;
                    }
				annealingcounter++;
				}
			else if(annealingcounter<maxoutputconfigs)
				{
				annealingcounter++; 
                for(unsigned int i=0;i<totalNumSpins;i++) configuration << spinArray[i] << " "; 
                configuration << energy << endl;
				}
			accratio=accepted_current/double(temproundcounter);
			temperatures << counter << '\t' << temperature << '\t' << accratio << endl;
			energies << energy << endl;
			accepted+=accepted_current; accepted_current=0;
			temproundcounter=0;
			temproundsteps=stepIncrease(temproundsteps);
			}
		
		now = chrono::high_resolution_clock::now();
		timeDiff = chrono::duration_cast<chrono::milliseconds>(now-begin).count();
		}
		
    // WRAP UP --------------------------------------------------------------------------------------------
    
    result << endl;
	result << "Simulation took " << timeDiff/60./1000 << " minutes" << endl;
	result << "Finished after " << counter << " steps" << endl;
	result << "Final temperature: " << temperature << endl;
	result << "Average acceptance ratio: " << accepted/double(counter) << endl;
	result << endl;
	result << "final energy: " << energy << endl;
	result << "best energy: " << bestenergy << endl;
    result << "final probability: " << energy*probabilityNormFactor << endl;
    result << "best probability: " << bestenergy*probabilityNormFactor << endl;
	result << endl;
    
    ofstream finconf("finconf.dat");
    saveConfig(spinArray.data(), finconf);
	ofstream bestconf("bestconf.dat");
    saveConfig(bestSpinArray.data(), bestconf);
    
    return 0;
}
