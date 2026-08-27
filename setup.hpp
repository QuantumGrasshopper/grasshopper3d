#include "utilities.hpp"

void initLoad(unsigned char grid[], int spinArray[]);
void initRandom(unsigned char grid[], int spinArray[], gsl_rng *RNG);
void saveConfig(const int spinArray[], const std::string& filename);
void initialize(unsigned char grid[], int spinArray[], gsl_rng *RNG, std::string initconf);
