/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the function that runs a CHARMM script */

// It must be included by CHARMM.h
#ifndef CHARMM_Loading_Status
#error CHARMM methods must be included by CHARMM.h
#endif

// Implement the method
void CHARMM::run_charmm_script (const string& script, const string& LABEL) {
    // The script should be written to this file
    string fileName1 = LABEL + "_input.inp";
    // And it's output should go here
    string fileName2 = LABEL + "_output.out";
    // Write the script to the file
    ofstream output; output.open(fileName1.c_str());
    // If the file failed to open, throw an error
    if (!output.is_open()) {
        string error = "Failed to open this CHARMM script for writing: "
                     + fileName1 + "\n";
        throw PANTZ_error (error);}
    // Write the script to the file
    output << script; output.close();
    // The command to run the calculations
    // string command = "/home/shared/rjp0029_lab/charmm/exec/gnu/charmm < "
    string command = string(CHARMM_exec) + " < "
                   + fileName1 + " > " + fileName2;
    // Run the calculations
    int i = system(command.c_str());
    // Throw an error if something went wrong
    if (i != 0) {
        string error = "CHARMM calculations failed\n";
        throw PANTZ_error (error);}
    // End the function
}
