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

public:
    RandomValues() : gen(rd()) {
        massDistribution = std::uniform_int_distribution<int>(0, 2);
        energyDistribution = std::uniform_int_distribution<int>(0, 29);
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
};

std::vector<double> createRandomMomenta(int numParticles) {
    std::vector<double> randomMomenta(numParticles);

    RandomValues randomValues;
    std::normal_distribution<double> distribution(0.5, 0.25);

    auto start = std::chrono::high_resolution_clock::now();
    std::generate(randomMomenta.begin(), randomMomenta.end(),
                  [&]() { return 1 + 4 * distribution(randomValues.gen); });
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Vectorized execution time: " << duration.count() << " ms" << std::endl;

    return randomMomenta;
}