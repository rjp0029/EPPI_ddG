/* Created by the Pantazes Lab at Auburn University.
 *
 * This file declares a Rosetta namespace for doing Rosetta calculations using
 * PANTZ software. It is meant to be included by Methods.h and it includes it's
 * own set of header files that implement the methods in the namespace. */

// Use a header guard to make sure this file is only included once
#ifndef Rosetta_Header_Guard
#define Rosetta_Header_Guard 1

// Make sure this file is being included by Methods.h
#ifndef Methods_Loading_Status
#error Rosetta.h must be included by Methods.h
#endif

// Declare the namespace
namespace Rosetta {

    // Like CHARMM calculations, Rosetta calculations will be done in whatever
    // folder the PANTZ outputs are being written to. The CHARMM file management
    // functions can handle those tasks and will be used here as needed.
    
    // A function to make the label for the calculations
    string make_file_label (vector<PROT::Protein>&);
    // Delete the files that were made during a Rosetta calculation
    void clean_up (const string&, const bool, const string);
    // Output a set of proteins for Rosetta calculations
    void proteins_for_rosetta(vector<PROT::Protein>&, const string&, 
                              const bool);
    // Make a Rosetta energy minimization flag file
    void minimization_flag_file (const string&);
    // movemap file for fixed residue minimization
    void create_movemap_file (const string&, const vector<int>&);
    // Load proteins after Rosetta calculations
    void proteins_from_rosetta (vector<PROT::Protein>&, const string&);
    // A Rosetta Energy Minimization
    void Energy_Minimization (vector<PROT::Protein>&, ofstream&, const string&);
    // fixed residue energy minimization
    void Energy_Minimization_fixed_res (vector<PROT::Protein>&, ofstream&, const string&, vector<int>&);
    // The Interface Analysis calculations
    string determine_interface (const string&, const vector<PROT::Protein>&);
    void interface_flag_file (const string&, const string&);
    void Interface_Analyzer (vector<PROT::Protein>&, ofstream&, const string&,
                             const string&);
    float dg_separated(const string&);
    // Generate the per-residue interaction energies
    void per_residue_flag_file (const string&);
    void Per_Residue (vector<PROT::Protein>&, ofstream&, const string&);
    // get the pose numbering
    int get_pose_numbering(string, char, int);
    int get_pose_numbering(PROT::PDB*, char, int);
    int get_pose_numbering(vector<PROT::Protein>&, char, int);
    // End the namespace
};

// Include the header files that implement those methods. Use a preprocessor
// command to make sure they're only included by this file
#define Rosetta_Loading_Status 1

// Include the Rosetta header files
#include "Rosetta/file_management.h"
#include "Rosetta/proteins.h"
#include "Rosetta/Minimization.h"
#include "Rosetta/Interface.h"
#include "Rosetta/Per_Residue.h"
#include "Rosetta/pose_numbering.h"

// Delete the loading status variable
#undef Rosetta_Loading_Status

// End the header guard from the start of the file
#endif
