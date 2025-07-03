/* Created by the Pantazes Lab at Auburn University.
 *
 * This is the main Protein Analysis ToolZ program. It is responsible for
 * reading in an instruction file and implementing the specified calculations it
 * requests. */

// Include the header files from the PANTZ source folder
#include "./source/Text.h"
#include "./source/PANTZ_error.h"
#include "./source/Proteins.h"
#include "./source/Interface.h"
#include "./source/Protocols.h"

// The program requires a command line argument
int main (int argc, char *argv[]) {
    // Put all the actual calculations in a try statement
    try {
        // The name of an instruction file must be provided in the command line
        if (argc != 2) {
            stringstream c1; c1 << argc-1;
            string error = "The PANTZ executable code requires exactly 1 command "
                           "line argument (the name of an instruction file), not "
                         + c1.str() + "\n";
            throw PANTZ_error(error);}
        
        // Get the name of the instruction file
        string instructionFile = argv[1];
        // Create the interface object that reads in that file and gets the
        // instructions
        Interface interface (instructionFile);

        // Validate that there are instructions (the interface should have done
        // this, but it is fine to repeat that check) and that the first
        // instruction specifies the type of calculation to run
        if (interface.commands() == 0) {
            string error = "No commands were found in " + instructionFile + "\n";
            throw PANTZ_error (error);}
        else if (interface.command_type(0) != "calculation type") {
            string error = "The first instruction in a PANTZ instruction file "
                           "must be a calculation type.\n";
            throw PANTZ_error (error);}

        // Get the specified calculation type
        string kind = (*interface.command(0))[0];
        // Convert it to lower case letters
        Text::lower(kind);
        // Run the appropriate program
        if (kind == "adhoc") {
            // Make the class object and run it's calculations
            PROTOCOL::ADHOC calculations (interface);
            calculations.run(interface);}
        // Otherwise, the calculation type is not supported and an error should
        // be thrown
        else {
            string error = "PANTZ does not support this calculation type: " 
                         + kind + "\n";
            throw PANTZ_error (error);}}
    
    // If any error occurred
    catch (PANTZ_error& e) {
        // Print the error message to the screen
        cout << e.what();
        // Return a value indicating that the program did not run correctly
        return 1;}

    // End the program
    return 0;
}
