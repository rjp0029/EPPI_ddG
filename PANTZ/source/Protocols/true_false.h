/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the true / false method for Protocols */

// This file is supposed to be included by Protocol.h
#ifndef Protocol_Loading_Status
#error General Protocol functions must be included by Protocol.h
#endif

// See if a string is an acceptable value to be interpreted as true or false
bool PROTOCOL::true_false (const string& value) {
    // Copy the value to a string to be edited
    string word = value;
    // Make the word lower case
    Text::lower(word);
    // If it is one of the favorable answers, change it to a standard
    if (((word == "yes") || (word == "y")) ||
        ((word == "true") || (word == "t"))) {word = "t";}
    // If it is unfavorable, convert it to a standard word
    else if (((word == "no") || (word == "n")) ||
             ((word == "false") || (word == "f"))) {word = "f";}
    // Otherwise, throw an error
    else {
        string error = value + " is not an acceptable value for a yes or no "
                       "instruction.\n";
        throw PANTZ_error (error);}
    // Return the appropriate value
    if (word == "t") {return true;}
    return false;
}
