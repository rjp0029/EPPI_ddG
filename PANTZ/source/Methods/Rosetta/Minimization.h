/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements energy minimization functions for Rosetta */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// Create a flag file for Rosetta calculations
void Rosetta::minimization_flag_file (const string& fileLabel) {
    // The contents are:
    string contents = "-in:file:s " + fileLabel + ".pdb \n"
                      "-run:min_type lbfgs_armijo_nonmonotone_atol\n"
                      "-run:min_tolerance 0.001\n"
                      "-out:path::all ./\n"
                      "-out:overwrite\n"
                      "-ignore_zero_occupancy false\n";
    // Write those contents to the flag file
    string fileName = fileLabel + "_flags.txt";
    ofstream output; output.open(fileName.c_str());
    if (!output.is_open()) {
        string error = "Failed to open " + fileName + " for Rosetta "
                       "calculations.\n";
        throw PANTZ_error (error);}
    output << contents;
    output.close();
    // End the function
}

// The function that actually runs a Rosetta energy minimization
void Rosetta::Energy_Minimization (vector<PROT::Protein>& proteins,
              ofstream& output, const string& path) {
    // If the output stream is open, include a message in it
    if (output.is_open()) {
        output << "Rosetta Energy Minimization started on "
               << METHODS::time_stamp() << endl;}
    // Get the current working directory
    string cwd = CHARMM::get_cwd(path);
    // Change to the output folder
    CHARMM::change_directory(path);
    // Make a label for naming everything
    string fileLabel = make_file_label (proteins);
    // Output the proteins
    proteins_for_rosetta(proteins, fileLabel, true);
    // Make the minimization flag file
    minimization_flag_file (fileLabel);
    // Create the command to run the minimization
    string command = string(ROSETTA_MIN_exec) + " @"
                   + fileLabel + "_flags.txt > " + fileLabel + "_output.out";
    // Run that command
    int i = system(command.c_str());
    // Load the proteins
    proteins_from_rosetta(proteins, fileLabel);
    // Clean up the created files
    clean_up (fileLabel, true, "minimization");
    // Move back to the original folder
    CHARMM::change_directory(cwd);
    // If appropriate, update the output file
    if (output.is_open()) {
        output << "Calculations ended on " + METHODS::time_stamp() 
               << "\n" << endl;}
}

void Rosetta::Energy_Minimization_fixed_res(vector<PROT::Protein>& proteins,
                                  ofstream& output, const string& path,
                                  vector<int>& fixed_residues) {
    // If the output stream is open, include a message in it
    if (output.is_open()) {
        output << "Rosetta Energy Minimization started on "
               << METHODS::time_stamp() << endl;
    }
    // Get the current working directory
    string cwd = CHARMM::get_cwd(path);
    // Change to the output folder
    CHARMM::change_directory(path);
    // Make a label for naming everything
    string fileLabel = make_file_label(proteins);
    // Output the proteins
    proteins_for_rosetta(proteins, fileLabel, true);
    // Make the minimization flag file
    minimization_flag_file(fileLabel);
    // Create the move map file
    create_movemap_file(fileLabel, fixed_residues);
    // Create the command to run the minimization
    string command = string(ROSETTA_MIN_exec) + " @" + fileLabel + "_flags.txt > " + fileLabel + "_output.out";
    // Run that command
    int i = system(command.c_str());
    // Load the proteins
    proteins_from_rosetta(proteins, fileLabel);
    // Clean up the created files
    clean_up(fileLabel, true, "minimization");
    // remove the move map file
    CHARMM::delete_file(fileLabel + "_movemap.txt");
    // Move back to the original folder
    CHARMM::change_directory(cwd);
    // If appropriate, update the output file
    if (output.is_open()) {
        output << "Calculations ended on " + METHODS::time_stamp()
               << "\n" << endl;
    }
}

// Helper function to create the move map file
void Rosetta::create_movemap_file(const string& fileLabel, const vector<int>& fixed_residues) {
    ofstream movemap_file(fileLabel + "_movemap.txt");
    // get the pose numbering
    if (movemap_file.is_open()) {
        // Set all residues to movable by default
        movemap_file << "RESIDUE * NO" << endl;
        // Set the specified residues to fixed
        for (int res : fixed_residues) {
            movemap_file << "RESIDUE " << res << " CHI" << endl;
        }
        movemap_file.close();
    }
    // Update the flags file to include the move map
    ofstream flags_file(fileLabel + "_flags.txt", ios_base::app);
    if (flags_file.is_open()) {
        flags_file << "-movemap " << fileLabel + "_movemap.txt" << endl;
        flags_file.close();
    }
}

