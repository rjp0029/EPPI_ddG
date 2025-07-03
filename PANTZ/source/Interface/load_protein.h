/* Created by the Pantazes Lab at Auburn University.
 *
 * This file processes multiple lines of commands to load a protein */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Implement the method
size_t Interface::load_protein (const size_t CI, vector<string>& line1) {
    // This flag is used to indicate if any error occurs
    bool errorFlag = true;
    // Confirm that the line1 vector has 2 parts and the 1st is to load a
    // protein. If this happens, a code-creator made an algorithm error. It is
    // not an issue with the command file, but with the code.
    if ((line1.size() != 2) || (line1[0] != "load protein")) {
        string error = "Algorithm error in Interface::load_protein.\n";
        throw PANTZ_error (error);}
    // There are several details that will be stored for this command. The first
    // one is a classification
    string classification = line1[1];
    // The others are the folder the PDB file is stored in. Default to the
    // current folder if none is provided
    string folder = "./";
    // The second is the file name
    string fileName = "";
    // The third is the name of the protein chain in the PDB file
    string NameInFile = "";
    // And the fourth is the name of the protein chain to be used in the
    // calculations
    string NameForUse = "";
    // The index of the last checked command
    size_t Index = CI;
    // Put this in a try statement so that errors for failing to find properly
    // formatted content can be caught
    try { 
        // Get the next command, split into parts
        vector<string> parts;
        Index = get_next_command(parts, Index+1, true);
        // If this is a folder specification, store that information. Providing
        // a folder is option and the current folder is the default if not
        // specified
        if (parts[0] == "folder") {
            // Store the folder
            folder = parts[1];
            // Clear the parts, then get the next command
            parts.clear(); Index = get_next_command(parts, Index+1, true);}
        // This should be a file specification
        if (parts[0] == "file") {
            fileName = parts[1];
            // Get the next command (chain)
            parts.clear(); Index = get_next_command(parts, Index+1, true);
            if (parts[0] == "chain") {
                // Get the string
                string how = parts[1];
                // Split it into pieces
                parts.clear(); Text::split(parts, how);
                // If there is 1 piece, it should be both protein names
                if (parts.size() == 1) {
                    // It should be a single character
                    if (parts[0].size() == 1) {
                        char L = parts[0][0];
                        // Validate it
                        CHECK::protein_name (L);
                        // Set the name as both
                        NameInFile = parts[0];
                        NameForUse = parts[0];
                        // All needed information has been found, so set the
                        // error flag to false
                        errorFlag = false;}}
                // if there are 3 pieces
                else if (parts.size() == 3) {
                    // Make the middle one lower case
                    Text::lower(parts[1]);
                    // It should be 'as'
                    if (parts[1] == "as") {
                        // Validate the 1st name
                        if (parts[0].size() == 1) {
                            char L = parts[0][0];
                            CHECK::protein_name(L);
                            NameInFile = parts[0];
                            // Validate the 2nd name
                            if (parts[2].size() == 1) {
                                L = parts[2][0];
                                CHECK::protein_name(L);
                                NameForUse = parts[2];
                                // All needed information has been found
                                errorFlag = false;}}}}}}
        // If the error flag is still set to true, there is a problem
        if (errorFlag) {
            string error = "Load Protein formatting error.\n";
            throw PANTZ_error (error);}}
    // If errors were triggered in the try statement, update them and re-throw
    // them.
    catch (PANTZ_error& e) {
        stringstream c1; c1 << CI + 1;
        string error = "This error occurred starting from the Load Protein "
                       "command in line " + c1.str() + "\n";
        throw PANTZ_error (e, error);}
    // To reach this point, no error occurred. Store the information
    vector<string> details;
    details.push_back(classification);
    details.push_back(folder);
    details.push_back(fileName);
    details.push_back(NameInFile);
    details.push_back(NameForUse);
    m_command_types.push_back(line1[0]);
    m_command_details.push_back(details);
    // Return the index of the last line that was evaluated
    return Index;
}
