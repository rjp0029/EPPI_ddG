/* Created by the Pantazes Lab at Auburn University.
 * 
 * This file implements a method that is a general strategy for processing an
 * instruction. */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// The general protocol for processing single line commands
void Interface::general_processing (const vector<string>& parts) {
    // store the command type
    m_command_types.push_back(parts[0]);
    // Make a vector that is just the command details
    vector<string> details; details.push_back(parts[1]);
    // Store it in the command details
    m_command_details.push_back(details);
}
