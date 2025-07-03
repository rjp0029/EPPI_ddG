/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the folder management functions that are needed for
 * working with CHARMM calculations. It would be nice to avoid these, but at the
 * moment there is no known way to the authors of this code to allow CHARMM to
 * not require only lower case letters in file names. */

// Make sure this file is being included from the CHARMM.h header file
#ifndef CHARMM_Loading_Status
#error CHARMM functions must be included from CHARMM.h
#endif

// This module is needed for changing directories
#include <unistd.h>

// Implement the delete file method
void CHARMM::delete_file (const string& fileName) {
    // Delete the file using a system command
    string command = "rm " + fileName;
    // Run the command
    int i = system(command.c_str());
    // Throw an error if it failed
    if (i != 0) {
        string error = "Failure to delete: " + fileName + "\n";
        throw PANTZ_error (error);}
}

// The overly convoluted method that was identified to get the string of the
// current working directory
string CHARMM::get_cwd (const string& outputPath) {
    // Create a unique file name for the content
    string fileName = outputPath + "PANTZ_CHARMM_temp_file.txt";
    // Create a system command to store the information that file
    string command = "pwd > " + fileName;
    // Execute that command
    int i = system(command.c_str());
    // If the command failed
    if (i != 0) {
        string error = "Failure to find the current working directory\n";
        throw PANTZ_error (error);}
    // Open the file with the information
    ifstream input; input.open(fileName.c_str());
    // If it failed, throw an error
    if (!input.is_open()) {
        string error = "Failure to find the current working directory\n";
        throw PANTZ_error (error);}
    // Get the line of information
    string cwd; getline(input, cwd); Text::strip(cwd); input.close();
    // Delete the file that stored it
    delete_file(fileName);
    // Return the current working directory
    return cwd;
}

// Change a directory during CHARMM calculations
void CHARMM::change_directory (const string& path) {
    // Use the unistd.h function to do so
    int i = chdir(path.c_str());
    // Throw an error if it failed
    if (i != 0) {
        string error = "Failure to change directory to: " + path + "\n";
        throw PANTZ_error (error);}
}

// Generate the name of a file for CHARMM calculations
string CHARMM::make_protein_name (PROT::Protein& prot, const bool input) {
    // Store the name here
    string name = "protein_";
    name += prot.name();
    if (input) {name += "_input.pdb";}
    else {name += "_output.pdb";}
    // Make sure everything is lower case (the protein's name wasn't)
    Text::lower(name);
    // Return that name
    return name;
}

// Delete all the files made for and by CHARMM calculations
void CHARMM::clean_up (const string& LABEL, vector<PROT::Protein>& prots,
                       const bool proteinsOutput, const bool keepScripts) {
    // Store the names of the files to delete in this vector
    vector<string> fileNames;
    // Deal with the input and output scripts
    string name = LABEL + "_input.inp"; fileNames.push_back(name);
    if (keepScripts) {
        string command = "cp " + name + " charmm_input.inp";
        int i = system(command.c_str());}
    name = LABEL + "_output.out"; fileNames.push_back(name);
    if (keepScripts) {
        string command = "cp " + name + " charmm_output.out";
        int i = system(command.c_str());}
    // The proteins' names
    for(size_t i=0; i<prots.size(); ++i) {
        // The name of the protein going into charmm
        name = make_protein_name(prots[i], true);
        fileNames.push_back(name);
        // If the protein was output from CHARMM
        if (proteinsOutput) {
            name = make_protein_name(prots[i], false);
            fileNames.push_back(name);}
    }
    // add the inp and out files
    fileNames.push_back("charmm_input.inp");
    fileNames.push_back("charmm_output.out");
    // fileNames.push_back("energy_minimization.inp");
    // fileNames.push_back("energy_minimization.out");
    // Delete all the files
    for(size_t i=0; i<fileNames.size(); ++i) {
        string command = "rm " + fileNames[i];
        int j = system(command.c_str());}
    // This function doesn't throw an error if the files aren't deleted
};
