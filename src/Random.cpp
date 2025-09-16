#include <ctime>
#include <random>
#include <vector>
#include <stdexcept>

#include "Random.h"

thread_local std::mt19937 rng_mt;
thread_local std::minstd_rand rng_minstd;
thread_local std::ranlux24 rng_ranlux24;

thread_local RNG_TYPE rng_type = MT;

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

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    switch(rng_type) 
    {
        case MT:
            return dist(rng_mt);

        case MINSTD:
            return dist(rng_minstd);
        
        case RANLUX24:
            return dist(rng_ranlux24);
        
        default:
            throw std::runtime_error("randomDouble01: Invalid RNG");
    }
}

double randomUniform01() {

    static thread_local std::mt19937 generator;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(generator);
}

double randomNormal(double mean, double stdDev) {

    static thread_local std::mt19937 generator;
    std::normal_distribution<double> dist(mean, stdDev);

    return dist(generator);
}

int randomInt(int low, int high) {

    if (high < low) {
        throw std::invalid_argument("randomInt: high is smaller than low");
    }

    std::uniform_int_distribution<int> dist(low, high);

    switch(rng_type)
    {
        case MT:
            return dist(rng_mt);
        
        case MINSTD:
            return dist(rng_minstd);
        
        case RANLUX24:
            return dist(rng_ranlux24);
        
        default:
            throw std::runtime_error("randomInt: Invalid RNG");
    }
}

double normalDist(double mean, double stdDev) {

    std::normal_distribution<double> dist(mean, stdDev);

    switch(rng_type)
    {
        case MT:
            return dist(rng_mt);
        
        case MINSTD:
            return dist(rng_minstd);
        
        case RANLUX24:
            return dist(rng_ranlux24);
        
        default:
            throw std::runtime_error("normalDist: Invalid RNG");
    }
}

std::vector<double> sample_truncated_gaussian_reject(double sigma, double R) {
    
    if (!(sigma > 0.0)) {
        throw std::invalid_argument("sample_truncated_reject: invalid sigma");
    }
    if (!(R >= 0.0)) {
        throw std::invalid_argument("sample_truncated_reject: invalid R");
    }

    const double R2 = R*R;

    while (true) {
        double x = normalDist(0.0, sigma);
        double y = normalDist(0.0, sigma);

        if (x*x + y*y <= R2) {
            return {x, y};
        }
    }
}