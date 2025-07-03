/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the declarations of the Protocols of the PANTZ software.
 * They are the main sets of calculations that PANTZ is routinley used to carry
 * out. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Protocols_Guard
#define Protocols_Guard 1

// Include the Proteins and Interface headers, which include other needed header
// files. Also include the Methods.
#include "Proteins.h"
#include "Interface.h"
#include "Methods.h"

// Declare a namespace for Protocol classes
namespace PROTOCOL {

    // The adhoc protocol is for a user-defined series of calculations
    class ADHOC;

    // Different Protocols need to do certain common things. Those functions are
    // defined in the PROTOCOL namespace and implemented in header files that
    // are included before the main Protocol header files are included.

    // Validate that a command has an appropriate amount of information
    void validate_command (vector<string> *, const size_t, const string&, 
                           const string&);
    // Determine if a value meets one of the true / false criteria
    bool true_false (const string&);
    // Fill the gaps in a Protein Object predicting the structure along the way
    void spasm (PROT::Protein *);
    // fill the gaps if the predicted structure is already known
    void align_and_splice_residues (PROT::Protein *, PROT::Protein *, bool terminal = false);
    // predict the ddG value for a mutation in an interface
    float EPPI_ddg (PROT::PDB*, string, string, string, bool);
    // run the EPPI_ddg protocol on an ensemble of structures (computed from backrub protocol rosetta)
    float EPPI_ddg_ensemble (PROT::PDB*, string, string, string, bool);

    // End the PROTOCOL namespace
};

// Use a pre-processor directive to make sure that protocols are all included
// from this file
#define Protocol_Loading_Status 1

// Include the common calculation header files
#include "Protocols/validate.h"
#include "Protocols/true_false.h"
#include "Protocols/spasm.h"
#include "Protocols/EPPI_ddg.h"

// Include the adhoc header file
#include "Protocols/ADHOC.h"

// Undefine the protocol loading status preprocessor variable
#undef Protocol_Loading_Status

// End the header guard from the start of the file
#endif
