/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements file management functions for Rosetta calculations */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// Make a label used in all files for a set of Rosetta calculations
string Rosetta::make_file_label (vector<PROT::Protein>& proteins) {
    // Start the label here
    string label = "rosetta_protein";
    // Make it plural when there are multiple proteins
    if (proteins.size() > 1) {label += "s";}
    // Add the lower case names of the proteins to the label
    for(size_t i=0; i<proteins.size(); ++i) {
        label += "_";
        label += proteins[i].name();}
    // Return the label
    return label;
}

// Delete the files made for and by Rosetta calculations
void Rosetta::clean_up (const string& label, const bool keepScripts, 
                        const string method) {
    // Store the names of the files in this vector
    vector<string> fileNames;
    // Store relevant names. First, the input PDB file
    string name = label + ".pdb"; fileNames.push_back(name);
    // There is always a flag file
    name = label + "_flags.txt"; fileNames.push_back(name);
    // There is always an output file
    name = label + "_output.out"; fileNames.push_back(name);
    // if keep files are being kept, copy that output to a different file
    if (keepScripts) {
        string command = "cp " + name + " rosetta_output.out";
        int i = system(command.c_str());}
    // if this is a minimization
    if (method == "minimization") {
        // There is an output proteins file
        name = label + "_0001.pdb"; fileNames.push_back(name);
        // There is a score file
        name = "score.sc"; fileNames.push_back(name);
        // if the score file should be kept
        if (keepScripts) {
            string command = "cp " + name + " rosetta_score.sc";
            int i = system(command.c_str());}}
    // If it is instead an interface analysis
    else if (method == "interface") {
        // There is a _RIA file
        name = label + "_RIA.out"; fileNames.push_back(name);
        // That file should always be kept
        string command = "cp " + name + " rosetta_interface_score.sc";
        int i = system(command.c_str());}
    // Or if it is a per-residue contribution calculation
    else if (method == "per residue") {
        // There is a _pr.out file
        name = label + "_pr.out"; fileNames.push_back(name);
        // That file should always be kept
        string command = "cp " + name + " rosetta_residue_scores.sc";
        int i = system(command.c_str());}
    // Otherwise, throw an error because this function doesn't know what to do
    else {
        string error = "Algorithm error: Rosetta::clean_up does not support "
                       "the " + method + " method.\n";
        throw PANTZ_error (error);}
    // Delete all the files
    for(size_t i=0; i<fileNames.size(); ++i) {
        string command = "rm " + fileNames[i];
        int j = system(command.c_str());}
    // This function doesn't throw an error if the files aren't deleted
};
