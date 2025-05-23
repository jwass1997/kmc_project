#pragma once 

enum RNG_TYPE {
    MT,
    MINSTD,
    RANLUX24
};

void setRandomSeed(int seed);

double randomDouble01();

int randomInt(int low, int high);

double normalDist(double mean, double stdDev);