/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the general open file methods for Protocols */

// This file is supposed to be included by Protocol.h
#ifndef Protocol_Loading_Status
#error General Protocol functions must be included by Protocol.h
#endif

// Implement the method for writing to files
void PROTOCOL::open_file (const string& fileName, ofstream& output, 
                          const bool overwrite) {
    // Make sure the output stream is not already open
    if (output.is_open()) {
        string error = "An already open output stream cannot be opened again.\n";
        throw PANTZ_error (error);}
    // If overwriting a previous file is not allowed
    if (!overwrite) {
        // Try to open a file with the specified name for reading
        ifstream input; input.open(fileName.c_str());
        // If it opened, throw an error
        if (input.is_open()) {
            input.close();
            string error = "Error opening a file for writing\n"
                         + fileName + " already exists and cannot be "
                           "overwritten.\n";
            throw PANTZ_error (error);}}
    // if the function reached this point, it is acceptable to open the file
    output.open(fileName.c_str());
    // If it didn't open, throw an error
    if (!output.is_open()) {
        string error = "Error opening a file for writing\n" + fileName
                     + " could not be opened.\n";
        throw PANTZ_error (error);}
}

// The method for reading from files
void PROTOCOL::open_file (const string& fileName, ifstream& input) {
    // If the input stream is already open, throw an error
    if (input.is_open()) {
        string error = "An already open input stream cannot be opened again.\n";
        throw PANTZ_error (error);}
    // Try to open the file
    input.open(fileName.c_str());
    // If the file didn't open, throw an error
    if (!input.is_open()) {
        string error = "Error opening a file for reading\n" + fileName
                     + " could not be opened.\n";
        throw PANTZ_error (error);}
}
