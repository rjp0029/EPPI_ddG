/* Created by the Pantazes Lab at Auburn University.
 *
 * This file processes commands for an Instruction file */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Implement the method
size_t Interface::process_command (size_t Index) {
    // Store parts of the command here
    vector<string> parts;
    // Find the next command. Don't throw an error if the end of the file is
    // reached
    Index = get_next_command (parts, Index, false);
    // If the end of the command file was reached, throw an error
    if (Index >= m_all_contents.size()) {return Index;}
    // Otherwise, evaluate the type of command that was found
    if (parts[0] == "calculation type") {general_processing(parts);}
    else if (parts[0] == "calculation name") {general_processing(parts);}
    else if (parts[0] == "output path") {general_processing(parts);}
    else if (parts[0] == "create output folder") {general_processing (parts);}
    else if (parts[0] == "create summary file") {general_processing(parts);}
    else if (parts[0] == "overwrite previous") {general_processing(parts);}
    else if (parts[0] == "load protein") {Index = load_protein(Index, parts);}
    else if (parts[0] == "write proteins") {general_processing(parts);}
    else if (parts[0] == "charmm") {general_processing(parts);}
    else if (parts[0] == "rosetta") {general_processing(parts);}
    // If the command is not recognized, throw an error
    else {
        string error = parts[0] + " is not a recognized instruction.\n";
        throw PANTZ_error (error);}
    // Return the Index of the last evaluated command
    return Index;
}
