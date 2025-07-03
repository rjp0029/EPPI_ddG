/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements interface analysis functions for Rosetta */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// The flag file needs to know how proteins interact in the interface. 
string Rosetta::determine_interface (const string& command,
                                     const vector<PROT::Protein>& proteins) {
    // Split the command into parts based on white space
    vector<string> parts; Text::split(parts, command);
    // There must be at least 2 parts
    if (parts.size() < 2) {
        string error = "Algorithm Error: invalid command for "
                       "Rosetta::determine_interface - " + command + "\n";
        throw PANTZ_error (error);}
    // Make sure the first 2 parts are interface and analysis
    Text::lower(parts[0]); Text::lower(parts[1]);
    if ((parts[0] != "interface") || (parts[1] != "analysis")) {
        string error = "Algorithm Error: invalid command given to the "
                       "Rosetta::Interface_Analyzer function - " + command 
                     + "\n";
        throw PANTZ_error (error);}
    // If there is only 1 protein, throw an error
    if (proteins.size() == 1) {
        string error = "The Rosetta Interface Analyzer function cannot run "
                       "for a system of only 1 protein.\n";
        throw PANTZ_error (error);}
    // Store the flag string here
    string flag = "";
    // If there were only 2 parts, there is a default behavior built into the
    // function
    if (parts.size() == 2) {
        // There must be exactly 2 proteins
        if (proteins.size() != 2) {
            stringstream c1; c1 << proteins.size();
            string error = "The proteins in an interface must be specified "
                           "for the Rosetta Interface Analyzer function for "
                           "systems of " + c1.str() + " proteins.\n";
            throw PANTZ_error(error);}
        // assemble the flag string
        flag += proteins[0].name();
        flag += "_";
        flag += proteins[1].name();
        return flag;}
    // These boolean flags will be used to indicate which proteins were
    // specified as part of the same group in an interface
    vector<bool> used; used.reserve(proteins.size());
    for(size_t i=0; i<proteins.size(); ++i) {used.push_back(false);}
    // Go through the parts of the command
    for(size_t i=2; i<parts.size(); ++i) {
        // If the size of the part is not 1, it can't be a protein name
        if (parts[i].size() != 1) {
            string error = "Invalid protein specifier in " + command + "\n";
            throw PANTZ_error (error);}
        // It must be upper case, so convert it
        Text::upper(parts[i]);
        // Get the character
        char L = parts[i][0];
        // Whether or not the protein named L is found
        bool found = false;
        // Find the protein
        for(size_t j=0; j<proteins.size(); ++j) {
            if (proteins[j].name() == L) {
                // If the protein was already used, throw an error
                if (used[j]) {
                    string error = "Duplicate protein specifier (";
                    error += L;
                    error += ") in " + command + "\n";
                    throw PANTZ_error (error);}
                // Mark the protein as used
                used[j] = true;
                // Mark the protein as found
                found = true;
                // Add the protein's name to the flag
                flag += L;
                // Stop the search
                break;}}
        // If the protein wasn't found, throw an error
        if (!found) {
            string error = "Protein ";
            error += L;
            error += " is not part of the system and cannot be part of "
                     "an analyzed interface - " + command + "\n";
            throw PANTZ_error (error);}
    }
    // Make sure there are proteins for the other part of the interface
    bool stillValid = false;
    for(size_t i=0; i<used.size(); ++i) {
        if (!used[i]) {stillValid = true; break;}}
    // If all proteins have been used, throw an error
    if (!stillValid) {
        string error = "Invalid command - all proteins were specified as part "
                       "of the same group for the interface analysis: "
                     + command + "\n";
        throw PANTZ_error (error);}
    // If the function reached this point, it can finish without issue
    flag += "_";
    for(size_t i=0; i<used.size(); ++i) {
        if (!used[i]) {
            flag += proteins[i].name();}}
    // REturn the flag
    return flag;
}

// Create a flag file for Rosetta calculations
void Rosetta::interface_flag_file (const string& fileLabel, 
                                   const string& interface) {
    // The contents are:
    string contents = "-in:file:s " + fileLabel + ".pdb \n"
                      "-interface " + interface + "\n"
                      "-compute_packstat=1\n"
                      "-add_regular_scores_to_scorefile=1\n"
                      "-out:file:score_only " + fileLabel + "_RIA.out\n"
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
void Rosetta::Interface_Analyzer (vector<PROT::Protein>& proteins,
              ofstream& output, const string& path, const string& command) {
    // If the output stream is open, include a message in it
    if (output.is_open()) {
        output << "Rosetta Interface Analyzer started on "
               << METHODS::time_stamp() << endl;}
    // Get the current working directory
    string cwd = CHARMM::get_cwd(path);
    // Change to the output folder
    CHARMM::change_directory(path);
    // Make a label for naming everything
    string fileLabel = make_file_label (proteins);
    // Make the interface flag
    string flag = determine_interface (command, proteins);
    // Output the proteins
    proteins_for_rosetta(proteins, fileLabel, true);
    // Make the flag file
    interface_flag_file (fileLabel, flag);
    // Create the command to run the calculations
    string what = string(ROSETTA_RIA_exec) + " @"
                + fileLabel + "_flags.txt > " + fileLabel + "_output.out";
    // Run that command
    int i = system(what.c_str());
    // Clean up the created files
    clean_up (fileLabel, true, "interface");
    // Move back to the original folder
    CHARMM::change_directory(cwd);
    // If appropriate, update the output file
    if (output.is_open()) {
        output << "Calculations ended on " + METHODS::time_stamp() 
               << "\n" << endl;}
}

float Rosetta::dg_separated(const string& file) {
    ifstream input; input.open(file.c_str());
    // split the second to last line by white space
    vector<string> parts;
    string line;
    while (getline(input, line)) {
        parts.clear();
        Text::split(parts, line);
    }
    // the 5th element is the binding energy
    float dg = atof(parts[5].c_str());
    return dg;
}
