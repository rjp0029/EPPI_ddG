/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is the header for functions and classes associated with using
 * Protein Data Bank - formatted information. It declares the PROT namespace,
 * and then includes other header files that implement the classes and methods
 * declared in that namespace. It also declares the CHECK namespace, which
 * contains functions to check that information meets PDB formatting
 * requirements. */

// Use a header guard to prevent this file from being included in a compiled
// program multiple times
#ifndef Proteins_Guard
#define Proteins_Guard 1

// Include other header files from the PANTZ directory
#include "Text.h"
#include "PANTZ_error.h"
#include "Macros.h"
// Include standard C++ files
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

// Declare the namespace for the classes
namespace PROT {
    // Define some constant information

    // PDB coordinates have 3 decimal places of significance, so a float is
    // appropriate to store that information. Use a typedef to allow that to
    // change later if needed
    typedef float coor;

    // The names of the 20 standard amino acids and 2 variants of histidine
    // that are commonly used
    const string AANames = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU "
                           "MET ASN PRO GLN ARG SER THR VAL TRP TYR "
                           "HSD HSE";
    // Split them into a vector
    const vector<string> AA3 = Text::split(AANames);
    // The corresponding 1 letter codes
    const string AACodes = "A C D E F G H I K L M N P Q R S T V W Y";
    const vector<string> AA1 = Text::split(AACodes);
    // The atoms that appear in the backbones of PDB-formatted amino acid
    // residues
    const string bbAtoms = "N H HN CA HA C O HN1 HN2 HT1 HT2 HT3 OT1 OT2";
    const vector<string> BackboneAtoms = Text::split(bbAtoms);
    // The number of coordinates an Atom has
    const size_t AtomCoordinates = 3;
    // The number of characters anticipated in a PDB Atom line
    const size_t AtomStringLength = 81;
    // The dielectric constant of water
    const float CCELEC = 331.843;

    // The matrix class is used to provide a standard container for linear
    // algebra-related tasks in this code. These are related to rotating and
    // moving proteins. While this could probably all be done with the Eigen
    // library's matrix class, it is not due to both historical,
    // code-development reasons and due to a desire to minimize reliance on
    // external code.
    class Matrix;
    
    // the KDtree class is used to store items in a way that
    // allows for fast searching of the items within a certain distance of
    // another or to get the k nearest neighbors to a point. This class is templated
    // to allow for different types of items to be stored in the tree. The items
    // must have a method to return the distance between two items.
    template<typename T>
    class KDtree;

    // The Atom class is a container of information about a single Atom in a PDB
    // file
    class Atom;
    // The AtomPtr class is a wrapper around a pointer to an Atom that makes it
    // easier to write some of the other expected code.
    //
    // This approach to classes should be fine, but does have the potential to
    // trigger a memory leak if improperly applied. That is something to keep an
    // eye on as PANTZ is developed.
    class AtomPtr;
    // Calculate the dihedral angle between 4 atoms. This function is
    // implemented at the end of the Atom.h header file
    coor calculate_dihedral (const Atom *, const Atom *,
                             const Atom *, const Atom *);
    coor calculate_dihedral (const Atom&, const Atom&, const Atom&, const Atom&);
    coor calculate_dihedral (AtomPtr&, AtomPtr&, AtomPtr&, AtomPtr&);

    // The Residue class is a container of Atoms that all are part of the same
    // amino acid. The ResiduePtr class is a wrapper that holds a pointer to a
    // Residue
    class Residue;
    class ResiduePtr;

    // A struct to store information about hydrogen bonds
    struct HydrogenBond;

    // a struct to store information about a salt bridge
    struct SaltBridge;

    // A struct to store information about a hydrophobic interaction
    struct Hydrophobic;

    // A Protein is all Residues that are part of the same primary structure
    class Protein;

    // PDB files can contain multiple copies of the same Protein. The Structure
    // class contains all of the copies of that same protein.
    class Structure;

    // The PDB class contains all of the information from a PDB-formatted file
    class PDB;

    // End the namespace
}

// A namespace for functions that validate that information follows PDB
// formatting specifications
namespace CHECK {
    // These functions will all throw errors if there is any problem
    void atom_number (const long);
    void atom_name (const string&);
    void alt_location (const char);
    void residue_name (const string&);
    void protein_name (const char);
    void residue_number (const long);
    void insertion_code (const char);
    void atom_coordinate (const PROT::coor);
    void occupancy (const float);
    void temperature (const float);
    void element (const string&);
    void charge (const string&);
    // Check to see whether or not a string is an amino acid name
    bool is_amino_acid (const string&);
    // Check to see whether or not a string is a backbone atom
    bool is_backbone_atom (const string&);
    // End the namespace
}

// Define a pre-processor variable. This is used to guarantee that the various
// header files are all included directly from this file.
#define Proteins_Loading_Status 1

// Include the header files that implement all of these classes
#include "PROT/Check.h"
#include "PROT/Matrix.h"
#include "PROT/KDtree.h"
#include "PROT/Atom.h"
#include "PROT/Residue.h"
#include "PROT/HydrogenBond.h"
#include "PROT/SaltBridge.h"
#include "PROT/Hydrophobic.h"
#include "PROT/Protein.h"
#include "PROT/Structure.h"
#include "PROT/PDB.h"

// Undefine the proteins loading status pre-processor variable
#undef Proteins_Loading_Status

// End the header guard from the start of the file
#endif
