#include <memory>
#include "TFile.h"
#include "TTree.h"

struct Data {
    double momentum;
    int typeOfParticle;
    double refractiveIndex;
    double pads[10][10];
};

void saveData() {
    // Create new file using TFile constructor
    auto file = std::make_unique<TFile>("data.root", "RECREATE");

    // Create a new TTree
    auto tree = std::make_unique<TTree>("T", "data tree");

    Data data;
    
    // Create branches in the tree
    tree->Branch("Momentum", &data.momentum, "momentum/D");
    tree->Branch("TypeOfParticle", &data.typeOfParticle, "typeOfParticle/I");
    tree->Branch("RefractiveIndex", &data.refractiveIndex, "refractiveIndex/D");
    tree->Branch("Pads", &data.pads, "pads[10][10]/D");

    // Fill tree with some data
    // In a real case, you would likely loop over some data collection here
    for (int i = 0; i < 100; i++) {
        data.momentum = i * 0.5;  // some dummy values
        data.typeOfParticle = i % 3;
        data.refractiveIndex = i * 0.1;
        for(int j=0; j<10; j++){
            for(int k=0; k<10; k++){
                data.pads[j][k] = i+j+k; // some dummy values
            }
        }
        tree->Fill();
    }

    // Write the tree into the file
    tree->Write();
}
