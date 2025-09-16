#pragma once
#include <random>
#include <vector>

enum RNG_TYPE {
    MT,
    MINSTD,
    RANLUX24
};

extern thread_local std::mt19937 rng_mt;
extern thread_local std::minstd_rand rng_minstd;
extern thread_local std::ranlux24 rng_ranlux24;

extern thread_local RNG_TYPE rng_type;

void setRandomSeed(int seed);

double randomDouble01();

double randomUniform01();

double randomNormal(double mean, double stdDev); 

int randomInt(int low, int high);

double normalDist(double mean, double stdDev);

std::vector<double> sample_truncated_gaussian_reject(double sigma, double R);