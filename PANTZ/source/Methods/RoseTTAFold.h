/* Created by the Pantazes Lab at Auburn University.
 *
 * This file declares a RoseTTAFold namespace. */

// Use a header guard to make sure this file is only included once
#ifndef RoseTTAFold_Header_Guard
#define RoseTTAFold_Header_Guard 1

// Make sure this file is being included by Methods.h
#ifndef Methods_Loading_Status
#error RoseTTAFold.h must be included by Methods.h
#endif

// Declare the namespace
namespace RoseTTAFold {

    // predict a structure using RoseTTAFold given a primary sequence and protein char
    PROT::Protein * predict_structure(string &, char &);

    // End the namespace
};

// Include the header files that implement those methods. Use a preprocessor
// command to make sure they're only included by this file
#define RoseTTAFold_Loading_Status 1

// Include the RoseTTAFold header files
#include "RoseTTAFold/predict_structure.h"

// Delete the loading status variable
#undef RoseTTAFold_Loading_Status

// End the header guard from the start of the file
#endif
