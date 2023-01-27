#include "common3d.hpp"


using namespace std;

unsigned int totalNumSpins = 5000;
double cellsize = pow(1/double(totalNumSpins),1./3.);
unsigned int gridSize = 50;
unsigned int gridVolume = gridSize*gridSize*gridSize;

int main() {
    
    // test parameters
    double d = 0.2;
    long unsigned int steps = 1000000;
    long unsigned int acceptance_counter=0;
    
    // RNG
	auto seed=12345;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
    
    bool grid[gridVolume];					//true if spin=1 at this grid point
    double energyGrid[gridVolume];
	int spinArray[totalNumSpins];				//grid point where any spin is
	int noSpinArray[gridVolume-totalNumSpins];		//complementary to above: grid point where no spin is
	vector< pair<int,double> > dNeighbourTemplate;
	vector< pair<int,double> > dNeighbourTable[gridVolume];	//for each grid point: list of grid points that are its d-neighbours with corresponding energies

    // construct generic neighbor list
    auto begin = chrono::high_resolution_clock::now();
    
    int center = gridSize*gridSize*gridSize/2+gridSize*gridSize/2+gridSize/2;
    double thisEnergyContribution;
    tuple<double,double,double> currentPosition=findPosition(center);
    for(unsigned int j=0;j<gridVolume;j++)
        {
        thisEnergyContribution=contributionEnergy(d,euclideanDistance(currentPosition,findPosition(j)));
        pair<int,double> thisPair(j,thisEnergyContribution);
        if(thisEnergyContribution > EPS) dNeighbourTemplate.push_back(thisPair);
        }
    
    // make list of neighbor lists for each grid cell (?)    
    //    int x = input%gridSize;
    //    int y = ((input - x)/gridSize)%gridSize;
    //    int z = (input - x - gridSize*y)/gridSize/gridSize;
	for(unsigned int j=0;j<dNeighbourTemplate.size();j++)
		{
        int coord = dNeighbourTemplate[j].first;
        int relx = xcoord(coord) - gridSize/2;
        int rely = ycoord(coord) - gridSize/2;
        int relz = zcoord(coord) - gridSize/2;
        double relativeEnergy = dNeighbourTemplate[j].second;
        for(unsigned int i=0;i<gridVolume;i++)
            {
            int gridLocationx = xcoord(i) + relx;
            int gridLocationy = ycoord(i) + rely;
            int gridLocationz = zcoord(i) + relz;
            if(gridLocationx >= 0 && gridLocationy >= 0 && gridLocationz >= 0 && gridLocationx < gridSize && gridLocationy < gridSize && gridLocationz < gridSize)
                {
                pair<int,double> thisPair(gridLocationx + gridLocationy*gridSize + gridLocationz*gridSize*gridSize, relativeEnergy);
                dNeighbourTable[i].push_back(thisPair);
                }
            }
		}
		
    // check that grid is right
    //for(unsigned int j=0;j<dNeighbourTable[405].size();j++) cout << dNeighbourTable[405][j].first << '\t' << dNeighbourTable[405][j].second << endl;
    
    auto now = chrono::high_resolution_clock::now();
	cout << "Time to construct neighbors list: " <<chrono::duration_cast<chrono::milliseconds>(now-begin).count() << "ms" << endl;

    // initialize grid with random configuration
    for(unsigned int i=0;i<gridVolume;i++)
        {
        grid[i]=false;
        }
    int newSpinCoord; unsigned int spincounter=0;
    while(spincounter<totalNumSpins)
        {
        bool create=true;
        while(create==true)
            {
            newSpinCoord=gsl_rng_uniform_int (RNG, gridVolume);
            create=grid[newSpinCoord];
            }
        grid[newSpinCoord]=true;
        spinArray[spincounter]=newSpinCoord;
        spincounter++;
        }
        
    unsigned int noSpinCounter=0;
	for(unsigned int i=0;i<gridVolume;i++)
		{
		if(grid[i]==false) {noSpinArray[noSpinCounter]=i; noSpinCounter++;}
        // Fill the energy grid
        energyGrid[i]=0;
        for(unsigned int j=0;j<dNeighbourTable[i].size();j++)
			{
			if(grid[dNeighbourTable[i][j].first]==true)
				{energyGrid[i] += dNeighbourTable[i][j].second;}
			}
		}

    // perform a set of X regular MC updates, let's say we only accept improvements
    for(unsigned int counter=0;counter<steps;counter++)
        {
        //select random spins to destroy and to create
		int destroy=gsl_rng_uniform_int (RNG, totalNumSpins);
		int oldSpinCoord=spinArray[destroy];
		int create=gsl_rng_uniform_int (RNG, gridVolume-totalNumSpins);
		int newSpinCoord=noSpinArray[create];
		
        //calculate energy difference
		double energyDifference=energyGrid[newSpinCoord]-energyGrid[oldSpinCoord];
        if(isAround(d,euclideanDistance(findPosition(newSpinCoord),findPosition(oldSpinCoord))))//NOTE more efficient to keep this check explicit
            {
            energyDifference -= contributionEnergy(d,euclideanDistance( findPosition(newSpinCoord),findPosition(oldSpinCoord) ));
            }
    
		if(energyDifference>=0)
			{
			grid[oldSpinCoord]=false; grid[newSpinCoord]=true;
			spinArray[destroy]=newSpinCoord; noSpinArray[create]=oldSpinCoord;
        
            for(unsigned int j=0;j<dNeighbourTable[oldSpinCoord].size();j++)
                {
                energyGrid[dNeighbourTable[oldSpinCoord][j].first] -= dNeighbourTable[oldSpinCoord][j].second;
                }
            for(unsigned int j=0;j<dNeighbourTable[newSpinCoord].size();j++)
                {
				energyGrid[dNeighbourTable[newSpinCoord][j].first] += dNeighbourTable[newSpinCoord][j].second;
                }
        
			acceptance_counter++;
			}
        }

    // benchmark the time  
    auto end = chrono::high_resolution_clock::now();
    cout << "Benchmark time: " << chrono::duration_cast<chrono::milliseconds>(end-now).count() << "ms" << endl;
    cout << "Acceptance ratio: " << acceptance_counter/double(steps) << endl;
    
    return 0;
}
