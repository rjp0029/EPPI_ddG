/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the Per Residue energy function for Rosetta */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// Make a flag file
void Rosetta::per_residue_flag_file (const string& fileLabel) {
    // The contents are:
    string contents = "-in:file:s " + fileLabel + ".pdb\n"
                      "-out:file:silent " + fileLabel + "_pr.out\n"
                      "-out:overwrite\n"
                      "-ignore_zero_occupancy false\n";
    // Write the contents to the flag file
    string fileName = fileLabel + "_flags.txt";
    ofstream output; output.open(fileName.c_str());
    if (!output.is_open()) {
        string error = "Failed to open " + fileName + " for Rosetta "
                       "calculations.\n";
        throw PANTZ_error (error);}
    output << contents;
    output.close();
}

// The function that runs the calculation
void Rosetta::Per_Residue (vector<PROT::Protein>& proteins,
                           ofstream& output, const string& path) {
    // If the output stream is open, include a message in it
    if (output.is_open()) {
        output << "Rosetta Per Residue Energy calculations started on "
               << METHODS::time_stamp() << endl;}
    // Get the current working directory
    string cwd = CHARMM::get_cwd(path);
    // Change to the output folder
    CHARMM::change_directory(path);
    // Make a label for naming everything
    string fileLabel = make_file_label (proteins);
    // Output the proteins
    proteins_for_rosetta(proteins, fileLabel, false);
    // Make the flag file
    per_residue_flag_file (fileLabel);
    // Create the command to run the calculation
    string command = string(ROSETTA_REB_exec) + " @"
                   + fileLabel + "_flags.txt > " + fileLabel + "_output.out";
    // Run that command
    int i = system(command.c_str());
    // Clean up the created files
    clean_up (fileLabel, true, "per residue");
    // Move back to the original folder
    CHARMM::change_directory(cwd);
    // If appropriate, update the output file
    if (output.is_open()) {
        output << "Calculations ended on " + METHODS::time_stamp() 
               << "\n" << endl;}
}
