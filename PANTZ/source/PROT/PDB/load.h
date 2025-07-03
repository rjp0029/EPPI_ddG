/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the load method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Load the contents of the PDB file into this file
void PROT::PDB::load () {
    // Assemble the name of the file 
    string fileName = m_folder + m_name;
    // Attempt to open that file
    ifstream input; input.open(fileName.c_str());
    // If the file is not open, throw an error
    if(!input.is_open()) {
        string error = "Failure to open:\nFile: " + m_name + "\nLocation: ";
        if ((m_folder == "./") || (m_folder == "")) {
            error += "Current Folder\n";}
        else {error += m_folder + "\n";}
        throw PANTZ_error (error);}
    // The lines will be stored here
    string line;
    // For NMR files, only lines containing the first model will be loaded
    bool model_flag = false;
    // Go through the file's lines
    getline(input, line);
    while(!input.eof()) {
        // Strip whitespace from the line
        Text::strip(line);
        // If the line starts with the word "model"
        if (Text::startswith(line, "MODEL")) {
            // Split the line into pieces
            vector<string> parts = Text::split(line);
            // If there are at least two pieces and the second piece is not
            // "1", set the flag to true. If it is 1, set the flag to false
            if (parts.size() > 1) {
                if (parts[1] == "1") {model_flag = false;}
                else {model_flag = true;}}}
        // If the line should be stored, do so
        if(!model_flag) {m_lines.push_back(line);}
        // Get the next line
        getline(input, line);}
    // Close the input file
    input.close();
    // If no contents were identified, raise an error
    if (m_lines.size() == 0) {
        string error = "No contents were identified in:\nFile: " + m_name
                     + "\nLocation: " + m_folder + "\n";
        throw PANTZ_error (error);}
    // End this function
}
