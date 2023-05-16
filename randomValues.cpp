#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>


class RandomValues {
private:

	static constexpr double nm2eV = 1239.842609;
	const double mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938; // masses in
	std::array<double, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};

	static constexpr double arrWaveLenDefault[30] = {
	  162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
	  182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
	  202, 204, 206, 208, 210, 212, 214, 216, 218, 220};

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
