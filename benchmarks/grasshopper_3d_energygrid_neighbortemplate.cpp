#include "common3d.hpp"

// Benchmarking code for the grasshopper
// Setup: square grid, default delta function discretization, random config 
//        fixed temperature=1, fixed number of MC updates
//        compare MC time after neighborslist initiated
//        also note time to create the neighborslist

// Storing energies of individual grid points: 3d version
// Not storing entire grid of neighbors list, but just the template 
// This is because for large grid sizes we run into memory issues
// The 2d version of this is grasshopper_neighbor_template.cpp, but here we also have the energy grid

using namespace std;

struct gridcell {
    int x;
    int y;
    int z;
    double energy;    
};

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
    
    bool grid[gridSize][gridSize][gridSize];					//true if spin=1 at this grid point
    double energyGrid[gridSize][gridSize][gridSize];
    
    vector< tuple<int,int,int> > spinArray(totalNumSpins);
	vector< tuple<int,int,int> > noSpinArray(gridVolume-totalNumSpins);
    
	vector< gridcell > dNeighbourTemplate;

    // construct generic neighbor list
    auto begin = chrono::high_resolution_clock::now();
    
    double thisEnergyContribution;    
    tuple<int,int,int> centerPosition=make_tuple(gridSize/2,gridSize/2,gridSize/2);        
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
                    this_cell.x = i-gridSize/2;
                    this_cell.y = j-gridSize/2;
                    this_cell.z = n-gridSize/2;
                    this_cell.energy = thisEnergyContribution;
                    dNeighbourTemplate.push_back(this_cell);
                    }
                }
            }
        }
    
    auto now = chrono::high_resolution_clock::now();
	cout << "Time to construct neighbors list: " <<chrono::duration_cast<chrono::milliseconds>(now-begin).count() << "ms" << endl;

    // initialize grid with random configuration
    for(unsigned int i=0;i<gridSize;i++) 
        {
        for(unsigned int j=0;j<gridSize;j++)
            {
            for(unsigned int n=0;n<gridSize;n++)
                {
                grid[i][j][n]=false;
                }
            }
        }
    int newSpinCoord; unsigned int spincounter=0;
    while(spincounter<totalNumSpins)
        {
        bool create=true;
        while(create==true)
            {
            newSpinCoord=gsl_rng_uniform_int (RNG, gridVolume);
            create=grid[xcoord(newSpinCoord)][ycoord(newSpinCoord)][zcoord(newSpinCoord)];
            }
        grid[xcoord(newSpinCoord)][ycoord(newSpinCoord)][zcoord(newSpinCoord)]=true;
        spinArray[spincounter]=make_tuple(xcoord(newSpinCoord),ycoord(newSpinCoord),zcoord(newSpinCoord));
        spincounter++;
        }
        
    unsigned int noSpinCounter=0;
    for(unsigned int i=0;i<gridSize;i++) 
        {
        for(unsigned int j=0;j<gridSize;j++)
            {
            for(unsigned int n=0;n<gridSize;n++)
                {
                if(grid[i][j][n]==false) {noSpinArray[noSpinCounter]=make_tuple(i,j,n); noSpinCounter++;}
                // Fill the energy grid
                energyGrid[i][j][n]=0;
                for(unsigned int k=0;k<dNeighbourTemplate.size();k++)
                    {
                    int x = dNeighbourTemplate[k].x  +i;
                    int y = dNeighbourTemplate[k].y  +j;
                    int z = dNeighbourTemplate[k].z  +n;
                    if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                        {
                        if(grid[x][y][z]==true)
                            {energyGrid[i][j][n]+=dNeighbourTemplate[k].energy;}
                        }    
                    }
                }
            }
        }
        
//             // check that this was done right by checking energy of initial system
//                         double energy=0;
//                         for(unsigned int i=0;i<totalNumSpins;i++)
//                             {
//                             for(unsigned int j=i;j<totalNumSpins;j++)
//                                 {
//                                 energy+=contributionEnergy(euclideanDistance(spinArray[i],spinArray[j]),d);
//                                 }
//                             }
//                         cout << "Full energy of initial system (counting all pairs of spins): " << energy << endl;
// 
//                         energy = 0;
//                         for(unsigned int i=0;i<gridSize;i++)
//                             for(unsigned int j=0;j<gridSize;j++)
//                                 for(unsigned int k=0;k<gridSize;k++)
//                                     {
//                                     energy += energyGrid[i][j][k]*grid[i][j][k];
//                                     }
//                         cout << "Initial energy using energy grid: " << energy/2 << endl;


    // perform a set of X regular MC updates, let's say we only accept improvements
    for(unsigned int counter=0;counter<steps;counter++)
        {
        //select random spins to destroy and to create
		int destroy=gsl_rng_uniform_int (RNG, totalNumSpins);
		tuple<int,int,int> oldSpinCoord=spinArray[destroy];
		int create=gsl_rng_uniform_int (RNG, gridVolume-totalNumSpins);
		tuple<int,int,int> newSpinCoord=noSpinArray[create];
		
        //calculate energy difference
        double energyDifference=energyGrid[get<0>(newSpinCoord)][get<1>(newSpinCoord)][get<2>(newSpinCoord)]-energyGrid[get<0>(oldSpinCoord)][get<1>(oldSpinCoord)][get<2>(oldSpinCoord)];
        if(isAround(d,euclideanDistance(newSpinCoord,oldSpinCoord)))//NOTE more efficient to keep this check explicit
            {
            energyDifference -= contributionEnergy(d,euclideanDistance( newSpinCoord,oldSpinCoord ));
            }
    
		if(energyDifference>=0)
			{
			grid[get<0>(oldSpinCoord)][get<1>(oldSpinCoord)][get<2>(oldSpinCoord)]=false; 
            grid[get<0>(newSpinCoord)][get<1>(newSpinCoord)][get<2>(newSpinCoord)]=true;
			spinArray[destroy]=newSpinCoord; noSpinArray[create]=oldSpinCoord;
        
            for(unsigned int j=0;j<dNeighbourTemplate.size();j++)
                {
                int x = dNeighbourTemplate[j].x  +get<0>(oldSpinCoord);
                int y = dNeighbourTemplate[j].y  +get<1>(oldSpinCoord);
                int z = dNeighbourTemplate[j].z  +get<2>(oldSpinCoord);
                if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                    {
                    energyGrid[x][y][z] -= dNeighbourTemplate[j].energy;
                    }
                x = dNeighbourTemplate[j].x  +get<0>(newSpinCoord);
                y = dNeighbourTemplate[j].y  +get<1>(newSpinCoord);
                z = dNeighbourTemplate[j].z  +get<2>(newSpinCoord);
                if(x >= 0 && y >= 0 && z >=0 && x < gridSize && y < gridSize && z < gridSize)
                    {
                    energyGrid[x][y][z] += dNeighbourTemplate[j].energy;
                    }
                }
        
			acceptance_counter++;
			}
        }
        
//                     // check that this was done right by checking energy of final system
//                         energy=0;
//                         for(unsigned int i=0;i<totalNumSpins;i++)
//                             {
//                             for(unsigned int j=i;j<totalNumSpins;j++)
//                                 {
//                                 energy+=contributionEnergy(euclideanDistance(spinArray[i],spinArray[j]),d);
//                                 }
//                             }
//                         cout << "Full energy of initial system (counting all pairs of spins): " << energy << endl;
// 
//                         energy = 0;
//                         for(unsigned int i=0;i<gridSize;i++)
//                             for(unsigned int j=0;j<gridSize;j++)
//                                 for(unsigned int k=0;k<gridSize;k++)
//                                     {
//                                     energy += energyGrid[i][j][k]*grid[i][j][k];
//                                     }
//                         cout << "Initial energy using energy grid: " << energy/2 << endl;
        

    // benchmark the time  
    auto end = chrono::high_resolution_clock::now();
    cout << "Benchmark time: " << chrono::duration_cast<chrono::milliseconds>(end-now).count() << "ms" << endl;
    cout << "Acceptance ratio: " << acceptance_counter/double(steps) << endl;
    
    return 0;
}
