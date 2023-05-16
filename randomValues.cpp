#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>

class RandomValues {
private:
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<int> massDistribution;
    std::uniform_int_distribution<int> energyDistribution;
    std::normal_distribution<double> momentumDistribution;
    
public:
    double momentum;
    double mass;
    double energy;
    float refractiveIndex;


    RandomValues() : gen(rd()) {
        massDistribution = std::uniform_int_distribution<int>(0, 2);
        energyDistribution = std::uniform_int_distribution<int>(0, 29);
        momentumDistribution = std::normal_distribution<double>(0.5, 0.25);
        
        momentum = getRandomMomentum();
        mass = getRandomMass();
        energy = getRandomEnergy();
        refractiveIndex = calculateRefractiveIndex();

    }

    double getRandomMomentum() {
        return 1 + 4 * momentumDistribution(gen);
    }

    double getRandomMass() {
        int index = massDistribution(gen);
        return masses[index];
    }

    double getRandomEnergy() {
        int index = energyDistribution(gen);
        double photonEnergy = nm2eV / arrWaveLenDefault[index];
        return photonEnergy;
    }

   float calculateRefractiveIndex() {
        float photonEnergy = energy;
        float k = 1.177 + (0.0172) * photonEnergy;
        return k;
    }
};

/*
example usage :
    int numObjects = 5;  // Number of objects to create
    std::vector<RandomValues> randomObjects(numObjects);

*/