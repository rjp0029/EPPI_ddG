/* Created by the Pantazes Lab at Auburn University.
 *
 * This function separates an instruction line into a command type and the
 * corresponding command details */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Implement the method
void Interface::create_parts (vector<string>& parts, const string& instruction,
                              const size_t I) {
    // Split the instruction on the distinguishing character
    Text::split(parts, instruction, ':');
    // This flag will indicate that there aren't any errors
    bool errorFlag = true;
    // If there are two parts, check them
    if (parts.size() == 2) {
        // Strip whitespace from the first part then make it lower case
        Text::strip(parts[0]); Text::lower(parts[0]);
        // If it isn't empty, strip whitespace from the second part but do NOT
        // change it's case
        if (parts[0].size() > 0) {
            Text::strip(parts[1]);
            // If it isn't empty, there isn't an error here
            if (parts[1].size() > 0) {errorFlag = false;}}}
    // If any form of error occurred
    if (errorFlag) {
        stringstream c1; c1 << I+1;
        string error = "Command formatting error in line " + c1.str() + "\n";
        error += m_all_contents[I];
        throw PANTZ_error (error);}
}
