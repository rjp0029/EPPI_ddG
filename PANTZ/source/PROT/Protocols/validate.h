/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the general command validation methods for Protocols */

// This file is supposed to be included by Protocol.h
#ifndef Protocol_Loading_Status
#error General Protocol functions must be included by Protocol.h
#endif

// Validate that a command meets expected criteria
void PROTOCOL::validate_command (vector<string> * command, const size_t N,
                                 const string& calcType, 
                                 const string& commandType) {
    // If the number of entries in the command doesn't match what is expected
    if (command->size() != N) {
        stringstream c1; c1 << N;
        stringstream c2; c2 << command->size();
        string error = calcType + " calculations expect " + commandType
                     + " instructions to have " + c1.str() + " components, "
                       "not " + c2.str() + "\n";
        throw PANTZ_error (error);}
}
