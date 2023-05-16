#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>

std::vector<double> createRandomMomenta(int numParticles) {
    std::vector<double> randomMomenta(numParticles);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.5, 0.25);

    auto start = std::chrono::high_resolution_clock::now();
    std::generate(randomMomenta.begin(), randomMomenta.end(),
                  [&]() { return 1 + 4 * distribution(gen); });
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Vectorized execution time: " << duration.count() << " ms" << std::endl;

    return randomMomenta;
}


static constexpr double arrWaveLenDefault[30] = {
    162, 164, 166, 168, 170, 172, 174, 176, 178, 180,
    182, 184, 186, 188, 190, 192, 194, 196, 198, 200,
    202, 204, 206, 208, 210, 212, 214, 216, 218, 220
};
static constexpr double nm2eV = 1239.842609;


static constexpr double photonEnDefault[30] = {
    0.130857 0.132233 0.133619 0.135014 0.136419 0.137834 0.139259 0.140693 0.142137 0.143591
    0.145054 0.146527 0.148009 0.149501 0.151003 0.152515 0.154036 0.155568 0.157109 0.15866
    0.160221 0.161792 0.163372 0.164962 0.166562 0.168171 0.16979 0.171419 0.173057 0.174705
};

double randomMass(std::mt19937& gen) {
    std::uniform_int_distribution<int> distribution(0, 2);
    int index = distribution(gen);
    return masses[index];
}

double randomEnergy(std::mt19937& gen) {
    std::uniform_int_distribution<int> distribution(0, 29);
    int index = distribution(gen);
    //double photonEnergy = nm2eV / arrWaveLenDefault[index];
    double photonEnergy = photonEnDefault(index);
    return photonEnergy;
}
