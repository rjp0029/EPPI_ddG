// this is the main file for the predict_ddg program from the PREPPI paper by Richard, A.C. et al. (2025)
// #include <iostream>
#include "PANTZ/source/Protocols.h"
// chrono
#include <chrono>

int main(int argc, char * argv[]) {
    if (argc != 5) {
        cout << "Usage: " << argv[0] << " <pdb_file> <mutation> <interface> <output_path>" << endl;
        return 1;
    }
    string pdb_file = argv[1];
    string mutation = argv[2];
    string interface = argv[3];
    string output = argv[4];
    // create the output directory
    system(("mkdir -p " + output).c_str());
    // Create a PDB object
    PROT::PDB pdb(pdb_file);

    // start a timer
    auto start = chrono::high_resolution_clock::now();

    // Predict the ddG value for a mutation in an interface
    float ddg = PROTOCOL::EPPI_ddg(&pdb, mutation, interface, output, false);

    // stop the timer
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "Time taken: " << duration.count() << " seconds" << endl;

    // print the predicted ddG value to console
    cout << "Predicted ddG value: " << ddg << endl;
    return 0;
}