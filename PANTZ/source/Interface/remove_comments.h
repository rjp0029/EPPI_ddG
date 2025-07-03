/* Created by the Pantazes Lab at Auburn University.
 *
 * This function removes comments from an instruction file command */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Remove comments from a command
string Interface::remove_comments (const string& input) {
    // Set the output as equal to that input
    string output = input;
    // Strip whitespace
    Text::strip(output);
    // If it is empty, return it now
    if (output.size() == 0) {return output;}
    // If it starts with the comment symbol, return an empty string
    if (output[0] == '#') {
        output = "";
        return output;}
    // Try to split the line based on that symbol
    vector<string> parts; Text::split(parts, output, '#');
    // Return the 1st part
    output = parts[0];
    return output;
}
