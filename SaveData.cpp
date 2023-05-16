#include <memory>
#include <vector>
#include "TFile.h"
#include "TTree.h"

struct DataInfo {
    double momentum;
    int typeOfParticle;
    double refractiveIndex;
    double pads[10][10];
};

class DataSaver {
private:
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTree> tree;
    DataInfo data;

public:
    DataSaver(const std::string& filename) {
        // Create new file using TFile constructor
        file = std::make_unique<TFile>(filename.c_str(), "RECREATE");

        // Create a new TTree
        tree = std::make_unique<TTree>("T", "data tree");

        // Create branches in the tree
        tree->Branch("Momentum", &data.momentum, "momentum/D");
        tree->Branch("TypeOfParticle", &data.typeOfParticle, "typeOfParticle/I");
        tree->Branch("RefractiveIndex", &data.refractiveIndex, "refractiveIndex/D");
        tree->Branch("Pads", &data.pads, "pads[10][10]/D");
    }

    void fillData(const std::vector<DataInfo>& dataVector) {
        // Fill tree with data from the vector
        for (const auto& item : dataVector) {
            data = item;
            tree->Fill();
        }
    }

    void save() {
        // Write the tree into the file
        tree->Write();
    }
};
