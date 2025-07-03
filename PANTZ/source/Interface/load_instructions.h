/* Created by the Pantazes Lab at Auburn University.
 *
 * This file loads the contents from an instructions file for a set of PANTZ
 * calculations */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Load the contents of the instructions file
void Interface::load_instructions () {
    // Open the file
    ifstream input; input.open(m_instruction_file.c_str());
    if (!input.is_open()) {
        string error = "Failure to open: " + m_instruction_file + "\n";
        throw PANTZ_error (error);}
    // Read in and store its lines
    string line; getline(input, line);
    while(!input.eof()) {
        m_all_contents.push_back(line);
        getline(input, line);}
    // Close the file
    input.close();
    // If it was empty, throw an error
    if (m_all_contents.size() == 0) {
        string error = m_instruction_file + " did not contain content.\n";
        throw PANTZ_error (error);}
}


