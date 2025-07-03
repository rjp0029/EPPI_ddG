/* Created by the Pantazes Lab at Auburn University.
 *
 * This file declares a CHARMM namespace for doing CHARMM calculations using
 * PANTZ software. It is meant to be included by Methods.h and it includes it's
 * own set of header files that implement the methods in the namespace. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef CHARMM_header_guard
#define CHARMM_header_guard 1

// Make sure this file is being included by Methods.h
#ifndef Methods_Loading_Status
#error CHARMM.h must be included by Methods.h
#endif

// Declare the namespace
namespace CHARMM {

    // CHARMM only works with lower-case file names and paths. To avoid that,
    // the CHARMM calculations are done in the folder of the PANTZ calculations.
    // As a result, some file and folder related functions exist. The first one
    // deletes a specific file
    void delete_file (const string&);
    // Get the current working directory. There is almost certainly a better way
    // to do this, but I couldn't figure it out in a small amount of time.
    string get_cwd (const string&);
    // Change the current working directory
    void change_directory (const string&);
    // Make the name of a protein's file
    string make_protein_name (PROT::Protein&, const bool);
    // Delete all of the files made for and from CHARMM. 
    void clean_up (const string&, vector<PROT::Protein>&, const bool, 
                   const bool);
    // Set the warning and bomb levels in a script
    void warning_bomb (string&);
    // Load topology and parameter files
    void load_topology_parameter (string&);
    // Prepare proteins / files for CHARMM use
    void proteins_for_charmm (string&, vector<PROT::Protein>&);
    // Text to add missing atoms
    void add_missing_atoms (string&);
    // Text to run an energy minimization
    void energy_minimization (string&, const bool, const bool);
    // Output proteins from a CHARMM script
    void output_proteins (string&, vector<PROT::Protein>&);
    // Load proteins from CHARMM calculations
    void proteins_from_charmm (vector<PROT::Protein>&);
    // The function that actually runs a charmm script
    void run_charmm_script (const string&, const string&);
    // A CHARMM energy minimization script
    void Energy_Minimization (vector<PROT::Protein>&, ofstream&, const string&,
                              const string&);
    // A CHARMM script just to add missing atoms
    void Missing_Atoms (vector<PROT::Protein>&, ofstream&, const string&);

    // End the namespace
};

// include the header files that implement those methods. Use a preprocessor
// variable to make sure they're only included by this file
#define CHARMM_Loading_Status 1

// Include the CHARMM header files
#include "CHARMM/file_management.h"
#include "CHARMM/proteins.h"
#include "CHARMM/script.h"
#include "CHARMM/run.h"
#include "CHARMM/Minimization.h"
#include "CHARMM/Missing_Atoms.h"

// Delete the loading status variable
#undef CHARMM_Loading_Status

// End the header guard from the start of the file
#endif
