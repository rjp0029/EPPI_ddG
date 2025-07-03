/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the general time stamp method*/

// This file is supposed to be included by Methods.h
#ifndef Methods_Loading_Status
#error General Method functions must be included by Methods.h
#endif

// This function is the only one that uses the ctime library
#include <ctime>

// Implement the method
string METHODS::time_stamp () {
    // Get the current date and time
    time_t present = time(0);
    // Convert that into a character array
    char * timeInfo = ctime(&present);
    // Convert that into a string, which the writer of this function (Dr. Robert
    // Pantazes) is more comfortable working with
    string info = timeInfo;
    // Split that into parts
    vector<string> parts; Text::split(parts, info);
    // Make the string
    string output = parts[0] + ", " + parts[1] + " " + parts[2] + ", "
                  + parts[4] + " at " + parts[3];
    return output;
}
