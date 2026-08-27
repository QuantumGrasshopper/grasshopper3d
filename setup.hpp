#include "utilities.hpp"

void initLoad(unsigned char grid[], int spinArray[]);
void initRandom(unsigned char grid[], int spinArray[], gsl_rng *RNG);
void saveConfig(int *spinArray, std::ofstream &filename);
void initialize(unsigned char grid[], int spinArray[], gsl_rng *RNG, std::string initconf);
