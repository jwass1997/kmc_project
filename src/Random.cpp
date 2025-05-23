#include <ctime>
#include <random>
#include <stdexcept>

#include "Random.h"

static std::mt19937 rng_mt;
static std::minstd_rand rng_minstd;
static std::ranlux24 rng_ranlux24;

static RNG_TYPE rng_type = MT;

bool setRngType(const RNG_TYPE rng_type_) {

    rng_type = rng_type_;

    return true;
}

void setRandomSeed(int seed) {

    switch(rng_type) 
    {
    case MT:
        rng_mt.seed(seed);
        break;

    case MINSTD:
        rng_minstd.seed(seed);
        break;
    
    case RANLUX24:
        rng_ranlux24.seed(seed);
        break;
    
    default:
        throw std::runtime_error("Invalid RNG");
    }   

}

double randomDouble01() {

    switch(rng_type) 
    {
        case MT:
            return std::generate_canonical<double, 32>(rng_mt);

        case MINSTD:
            return std::generate_canonical<double, 32>(rng_minstd);
        
        case RANLUX24:
            return std::generate_canonical<double, 32>(rng_ranlux24);
        
        default:
            throw std::runtime_error("randomDouble01: Invalid RNG");
    }
}

int randomInt(int low, int high) {

    if (high < low) {
        throw std::invalid_argument("randomInt: high is smaller than low");
    }

    switch(rng_type)
    {
        case MT:
            return std::uniform_int_distribution<int>(low, high)(rng_mt);
        
        case MINSTD:
            return std::uniform_int_distribution<int>(low, high)(rng_minstd);
        
        case RANLUX24:
            return std::uniform_int_distribution<int>(low, high)(rng_ranlux24);
        
        default:
            throw std::runtime_error("randomInt: Invalid RNG");
    }
}

double normalDist(double mean, double stdDev) {

    switch(rng_type)
    {
        case MT:
            return std::normal_distribution(mean, stdDev)(rng_mt);
        
        case MINSTD:
            return std::normal_distribution(mean, stdDev)(rng_minstd);
        
        case RANLUX24:
            return std::normal_distribution(mean, stdDev)(rng_ranlux24);
        
        default:
            throw std::runtime_error("normalDist: Invalid RNG");
    }
}