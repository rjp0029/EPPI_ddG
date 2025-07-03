/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains an Interface class that is intended to do standard
 * processing of script files for PANTZ calculations. */ 

// Use a header guard to make sure the class is only included in a compiled
// program a single time
#ifndef Interface_Check
#define Interface_Check 1

// Include the Text and PANTZ_error header files
#include "Text.h"
#include "PANTZ_error.h"
// Include standard C++ libraries for input/output purposes
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

// Declare the Interface class
class Interface {

    // The information stored in the class is private
    private:
        // The name of the file containing the instructions
        string m_instruction_file;
        // The contents of that file
        vector<string> m_all_contents;
        // The things the program is supposed to do. These are a classification
        // and then a processed set of what to do as created by this class
        vector<string> m_command_types;
        vector<vector<string> > m_command_details;

    // Private methods that control class behavior
    private:
        // Load the instructions from a file
        void load_instructions ();
        // Remove comments from a command
        string remove_comments (const string&);
        // Separate it into parts
        void create_parts (vector<string>&, const string&, const size_t);
        // Find the next command line in the instructions
        size_t get_next_command (vector<string>&, size_t, const bool);
        // A general procedure for processing commands
        void general_processing (const vector<string>&);
        // How to process a command to load a protein
        size_t load_protein (const size_t, vector<string>&);
        // Process a specified command in the file
        size_t process_command (size_t);

    // The public methods of the Interface class
    public:
        // The class constructor
        Interface (const string&);
        // Information about the contents of that file
        size_t lines () const {return m_all_contents.size();}
        string line (const size_t) const;
        // Information about commands
        size_t commands () const {return m_command_types.size();}
        string command_type (const size_t) const;
        vector<string> * command (const size_t);
        // Write the full instruction file to a new file
        void write (const string&) const;

    // End the class definition
};

// Load the methods of the Interface class
#define PANTZ_Interface_Loading_Status 1

#include "Interface/load_instructions.h"
#include "Interface/remove_comments.h"
#include "Interface/create_parts.h"
#include "Interface/get_next_command.h"
#include "Interface/general_processing.h"
#include "Interface/load_protein.h"
#include "Interface/process_command.h"

// Undefine the preprocessor command that is used to make sure the files are
// being loaded from this file
#undef PANTZ_Interface_Loading_Status

// The public class methods

// The class constructor
Interface::Interface (const string& fileName) {
    // Store the instruction file name
    m_instruction_file = fileName;
    // Load the instructions
    load_instructions();
    // Process the commands
    size_t Index = 0;
    while (Index < m_all_contents.size()) {
        Index = process_command(Index) + 1;}
}

// Access to a specific line, instruction, or command
string Interface::line (const size_t i) const {
    if (i >= m_all_contents.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_all_contents.size();
        string error = c1.str() + " is not a valid line index in an instruction"
                       " file of " + c2.str() + " lines.\n";
        throw PANTZ_error (error);}
    return m_all_contents[i];
}

// Access to a command type
string Interface::command_type (const size_t i) const {
    if (i >= m_command_types.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_command_types.size();
        string error = c1.str() + " is not a valid command index in a file of "
                     + c2.str() + " commands.\n";
        throw PANTZ_error (error);}
    return m_command_types[i];
}

// Access to the details of a specific command
vector<string> * Interface::command (const size_t i) {
    if (i >= m_command_details.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_command_details.size();
        string error = c1.str() + " is not a valid command index in a file of "
                     + c2.str() + " commands.\n";
        throw PANTZ_error (error);}
    return &(m_command_details[i]);
}

// Write the instructions to a file
void Interface::write (const string& outputName) const {
    // Open the output file
    ofstream output (outputName.c_str());
    if (!output.is_open()) {
        string error = "Failed to open: " + outputName + "\n";
        throw PANTZ_error (error);}
    // Write the contents of the instruction file
    for(size_t i=0; i<m_all_contents.size(); ++i) {
        output << m_all_contents[i];}
    output.close();
}

// End the header guard from the start of the file
#endif
