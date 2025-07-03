/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the general make folder method */

// This file is supposed to be included by Methods.h
#ifndef Methods_Loading_Status
#error Methods::make_folder.h must be included by Methods.h
#endif

// Implement the method
void METHODS::make_folder (const string& folder, const bool overwrite) {
    // The system command to make the folder is
    string command = "mkdir " + folder;
    // Use the system function to run that command
    int i = system(command.c_str());
    // If i is 0, the command worked and the folder was made. Otherwise, it did
    // not run correctly. If the value is not 0 and overwriting is not allowed,
    // throw an error
    if ((i != 0) && (!overwrite)) {
        string error = "The program failed to make this folder: " + folder
                     + "\nIt likely already exists.\n";
        throw PANTZ_error (error);}
}                   
