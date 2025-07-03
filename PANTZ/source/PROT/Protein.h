/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the definitions of the Protein class, which is a container
 * of Residues, for working with PDB-formatted data. It also includes all of the
 * header files that implement the class methods. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_Protein_Guard
#define Proteins_Protein_Guard 1

// Make sure the file is being included from Proteins.h
#ifndef Proteins_Loading_Status
#error Protein.h must be included by Proteins.h
#endif

// Confirm that the Check, Matrix, Atom, and Residue headers have been included.
// Because Residue does that for Check, Matrix and Atom, only Residue needs to
// be checked here
#ifndef Proteins_Residue_Guard
#error Protein.h must be included after Residue.h
#endif

// Declare the Protein class
class PROT::Protein {

    // The PDB class is a friend
    friend class PROT::PDB;

    // The information stored in the Protein is private
    private:
        // The Protein's name
        char m_name;
        // The number of residues
        size_t m_count;
        // The Residues themselves
        Residue * m_residues;

    // Private methods that control class behavior
    private:
        // Assign default values to class methods
        void initialize ();
        // Delete dynamically allocated memory
        void clean_up ();
        // Copy information from another instance of the class
        void copy (const Protein *);
        void copy (const Protein& other) {copy(&other);}
        // Validate that a set of Residues are acceptable for use in a Protein
        void check_residues (vector<Residue>&) const;
        // Private methods for moving and rotating a protein
        void private_move(const Matrix *, const char);
        void private_rotate(const Matrix *);

    // The public interface of the class
    public:
        // The class destructor
        ~Protein () {clean_up();}
        // The default constructor
        Protein () {initialize();}
        // A protein can be loaded from a vector of residues
        void load (vector<Residue>&);
        // Or a vector of atoms
        void load (vector<Atom>&);
        // Or a vector of strings
        void load (const vector<string>&);
        // Or a file
        void load (const string&, const string);
        // A Protein can also be initialized from those same inputs
        Protein (vector<Residue>& residues) {initialize(); load(residues);}
        Protein (vector<Atom>& atoms) {initialize(); load(atoms);}
        Protein (const vector<string>& contents) {initialize(); load(contents);}
        Protein (const string&, const string);
        // Copy construction of a Protein
        Protein (const Protein& other) {initialize(); copy(other);}
        Protein (const Protein * other) {initialize(); copy(other);}
        // Copy assignment
        void operator= (const Protein& other) {copy(other);}
        void operator= (const Protein * other) {copy(other);}
        // A formatted string of text containing all of the Protein's
        // information
        string str (const bool);
        // The number of atoms in the protein
        size_t number_of_atoms () const;
        // The number of the last atom in the protein
        long last_atom_number () const;
        // The INTERNAL number of the last residue in the protein. Internal
        // numbering starts at 1 and increases sequentially, with continuation
        // between proteins
        long last_residue_number () const;
        // The number of residues in the protein
        size_t size () const {return m_count;}
        // The name of the protein
        char name () const {return m_name;}
        // Change the internal numbering of the residues
        long renumber_residues (long);
        // Renumber the atoms
        long renumber_atoms (long);
        // Change the protein's name
        void set_name (const char);
        // Access to a Residue in the Protein using detailed name information
        Residue * operator () (const long, const char, const bool) const;
        // Note that it is not permitted to overload the () operator to have the
        // same command list and return a ResiduePtr instead. Which is OK, as
        // the ResiduePtr class should likely be used when nesting brackets
        // (e.g. protein[i][j][k]) 
        // Access to a ResiduePtr by index in this protein
        ResiduePtr operator[] (const size_t);
        // Move the Protein
        void move (const Matrix *, const char);
        void move (const Matrix *, const bool);
        void move (const Matrix&, const char);
        void move (const Matrix&, const bool);
        // Rotate the protein
        void rotate (const Matrix *);
        void rotate (const Matrix&);
        // Get a list of Atoms from the Protein
        void select_atoms (vector<Atom *>&, const string);
        void select_atoms (vector<AtomPtr>&, const string);
        // Center a Protein so the indicated atoms are at the origin on average
        Matrix center (const string);
        // Calculate a completeness score for a protein. This value is between 0
        // and 1 and is used in the context of PDB files
        double score () const;
        // Adjust CHARMM histidine names
        void for_charmm_histidine_fix ();
        void from_charmm_histidine_fix ();
        // Create a string of formatted text for the Protein for use in Rosetta
        string rosetta_str (long&, long&);
        // Fix atom naming after Rosetta
        void update_atoms_after_Rosetta ();
        // Get the FASTA sequence of the Protein
        string fasta (const size_t) const;
        // Calculate a Protein's dihedral angles
        void calculate_dihedrals ();
        // Make a duplicate of the Protein
        Protein duplicate () const;
    // End the class definition
};

// Define a loading status variable to make sure Protein methods are only
// included while this file is being loaded
#define ProteinClass_Loading_Status 1

// Include methods of the Protein class
#include "Protein/initialize.h"
#include "Protein/clean_up.h"
#include "Protein/copy.h"
#include "Protein/check_residues.h"
#include "Protein/move.h"
#include "Protein/rotate.h"
#include "Protein/load.h"
#include "Protein/renumber.h"
#include "Protein/str.h"
#include "Protein/numbers.h"
#include "Protein/operators.h"
#include "Protein/name.h"
#include "Protein/select_atoms.h"
#include "Protein/center.h"
#include "Protein/score.h"
#include "Protein/histidine.h"
#include "Protein/rosetta.h"
#include "Protein/fasta.h"
#include "Protein/dihedrals.h"
#include "Protein/duplicate.h"

// Undefine the Protein class loading guard
#undef ProteinClass_Loading_Status

// End the header guard from the start of the file
#endif
