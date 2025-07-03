/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains declarations of Methods of the PANTZ software. These are
 * common actions codes might take, such as loading a Protein or running an
 * energy minimization. More complex actions are Protocols. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef PANTZ_Methods_Guard
#define PANTZ_Methods_Guard 1

// Include the Proteins header file, which includes several other header files
// of interest
#include "Proteins.h"

// Methods are grouped into several categories, based on broad commonalities
// they have. Each category has a separate namespace

// Declare a namespace for general methods
namespace METHODS {

    // Create a string with the current date and time
    string time_stamp ();
    // A function to load a protein
    void load_protein (vector<string> *, vector<PROT::Protein>&, 
                       vector<PROT::PDB>&, ofstream&);
    // A function to output a protein
    void write_proteins (vector<string> *, vector<PROT::Protein>&,
                         const string&, ofstream&);
    // Make a folder to store files in
    void make_folder (const string&, const bool);
    // Open a file for writing
    void open_file (const string&, ofstream&, const bool);
    // Open a file for reading
    void open_file (const string&, ifstream&);
    // list directories
    vector<string> listdir(string);
    
    // A function to align the global structure of two proteins
    void global_align (PROT::Protein*, PROT::Protein*, const string&);
    // A function to align the local structure of two proteins specified by a set of atoms
    void align (PROT::Protein*, PROT::Protein*, vector<vector<PROT::Atom*>>);
    // set the atom sasa points in a given kdtree of atoms
    void set_sasa_points(vector<PROT::Atom*>);
    // method to make a mutation
    PROT::PDB make_mutation(PROT::PDB*, string&, string&);
    vector<PROT::Protein> make_mutation(vector<PROT::Protein>&, string&, string&);
};

// Use a pre-processor directive to make sure that Methods are all included
// from this file
#define Methods_Loading_Status 1

// Include the relevant header files
#include "Methods/time_stamp.h"
#include "Methods/load_protein.h"
#include "Methods/write_proteins.h"
#include "Methods/make_folder.h"
#include "Methods/open_file.h"
#include "Methods/listdir.h"
#include "Methods/align.h"
#include "Methods/CHARMM.h"
#include "Methods/Rosetta.h"
#include "Methods/RoseTTAFold.h"
#include "Methods/set_sasa_points.h"
#include "Methods/mutation.h"
#include "Methods/EPPI.h"

// Undefine the protocol loading status preprocessor variable
#undef Methods_Loading_Status

// End the header guard from the start of the file
#endif
