/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the definition of the Residue class, which is a container
 * of Atoms, for working with PDB-formatted data. It also includes all of the
 * header files that implement the methods of the class */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_Residue_Guard
#define Proteins_Residue_Guard 1

// Make sure that this file is being included by the Proteins header file
#ifndef Proteins_Loading_Status
#error Residue.h must be included by Proteins.h
#endif

// Confirm that the Check, Matrix and Atom header files are all included.
// Because Atom does that confirmation for Check and Matrix, Residue only needs
// to check that Atom has been loaded
#ifndef Proteins_Atom_Guard
#error Residue.h must be included after Atom.h
#endif

// Define the Residue class
class PROT::Residue {

    // The Protein class is a friend
    friend class PROT::Protein;
    // The PDB class is also a friend
    friend class PROT::PDB;

    // The information stored in the class is private
    private:
        // The Atoms in the Residue
        PROT::Atom * m_atoms;
        // The number of them
        size_t m_count;
        // The residue's name
        string m_name;
        // The residue's number
        long m_number;
        // The residue's internal number in a protein or list of residues
        long m_internal;
        // The residue's insertion code (used when two distinct residues have
        // the same number)
        char m_insertion;
        // The residue's protein
        char m_protein;
        // Whether or not the Residue is the N-terminus of a protein
        bool m_N_terminus = false;
        // Whether or not the Residue is the C-terminus of a protein
        bool m_C_terminus = false;

        // The following attributes are only meaningful in the context of PDB
        // files
        // Whether or not the Residue is present in the PDB file
        bool m_present;
        // The number of Atoms in the Residue that are absent from a PDB file
        size_t m_missing_atoms;

        // Dihedral angles for amino acids in proteins
        coor m_phi;
        coor m_psi;
        coor m_omega;

        // the stbaility of the residue in the context of the protein
        bool m_prestable = false;
        // the stability of the residue in the context of the protein complex
        bool m_bound_stable = false;
        // the rotamers of the residue
        vector<PROT::Residue> m_rotamers;
        // the intraprotein neighbors of the residue
        vector<PROT::Residue*> m_intra_neighbors;
        // the interprotein neighbors of the residue
        vector<PROT::Residue*> m_inter_neighbors;
        
    // Private functions that control the behaviour of the Residue
    private:
        // Assign default values to the Residue's attributes
        void initialize ();
        // Clean up dynamically allocated memory
        void clean_up ();
        // Copy information from another Residue
        void copy (const Residue *);
        void copy (const Residue& other) {copy(&other);}
        // Private functions to set residue information
        void private_set_number (const long, const char, const bool);
        void private_set_protein (const char);
        // Check that a list of Atoms are acceptable for use in a Residue
        void check_atoms (const vector<PROT::Atom *>&) const;
        // The private rotate and move functions of a Residue that do not
        // error check the matrix
        void private_move (const Matrix *, const char);
        void private_rotate (const Matrix *);
        // size_t free_rotamers(vector<PROT::Residue>, vector<PROT::Residue*>);
        // the utility function to determine the stability of the residue
        // bool stable (vector<PROT::Residue*>, vector<PROT::Residue*>);

    // The public interface of the Residue class
    public:
        // The class destructor
        ~Residue () {clean_up();}
        // The default class constructor
        Residue () {initialize();}
        // Update or entirely create the Residue from a vector of Atoms
        void load (const vector<PROT::Atom *>&, const bool, const bool);
        void load (vector<AtomPtr>&, const bool, const bool);
        // Construct a Residue from a vector of Atoms
        Residue (const vector<PROT::Atom *>&);
        // Construct a residue from a vector of AtomPtrs
        Residue (vector<AtomPtr>&);
        // A Residue can be constructed from a string (its name) and a
        // character (the protein's name). This should only be used in PDB
        // files.
        Residue (const string&, const char);
        // Copy construction of a Residue
        Residue (const Residue& other) {initialize(); copy(&other);}
        Residue (const Residue * other) {initialize(); copy(other);}
        // Copy assignment
        void operator= (const Residue& other) {copy(&other);}
        void operator= (const Residue * other) {copy(other);}
        // Access to the Residue's information
        size_t size () const {return m_count;}
        string name () const {return m_name;}
        char AA1 () const;
        long number () const {return m_number;}
        long internal_number () const {return m_internal;}
        char insertion_code () const {return m_insertion;}
        char protein () const {return m_protein;}
        bool is_present () const {return m_present;}
        size_t missing_atoms () const {return m_missing_atoms;}
        coor phi () const;
        coor psi () const;
        coor omega () const;
        // Access to the Residue's Atoms by number or name
        PROT::Atom * get_atom (const size_t);
        PROT::Atom * get_atom (const string&);
        PROT::AtomPtr operator[] (const size_t);
        PROT::AtomPtr operator[] (const string&);
        // Functions that allow for the setting of a Residue's number and
        // protein information
        void set_number (const long, const char, const bool);
        void set_protein (const char);
        void rename (const string&);
        // Renumber the atoms in the Residue
        long renumber_atoms (long);
        // The number of the last Atom in the Residue
        long last_atom_number () const;
        // A string of all of the Residue's information
        string str (const bool);
        // Move the Residue, with error checking of the provided matrix
        void move (const Matrix *, const char);
        void move (const Matrix *, const bool);
        void move (const Matrix&, const char);
        void move (const Matrix&, const bool);
        void rotate (const Matrix *);
        void rotate (const Matrix&);
        // Select a subset of Atoms from the Residue
        void select_atoms (vector<Atom *>&, const string);
        void select_atoms (vector<AtomPtr>&, const string);
        // Move a Residue so that its center of mass is at the origin
        PROT::Matrix center (const string);
        // Position a Residue for rotamer packaging
        void position (PROT::Matrix&, const bool);

        // Calculate a "completeness" score between 0 and 1 for the Residue - it
        // is the fraction of atoms that are present that should be present
        double score () const;
        // Whether or not the Residue is an amino acid
        bool is_amino_acid () const {return CHECK::is_amino_acid(m_name);}

        // A function that corrects histidine residue names for use in CHARMM
        void for_charmm_histidine_fix ();
        // And a function that moves them back to HIS after CHARMM
        void from_charmm_histidine_fix ();

        // Create a string of text that represents the non-hydrogen atoms of the
        // residue for use in Rosetta
        string rosetta_str (long&, long&, const bool);
        // Update atom names after Rosetta calculations
        void update_atoms_after_Rosetta (const bool);
        // This function creates a duplicated copy of the Residue and all of
        // its atoms.
        PROT::Residue duplicate () const;

        // x y and z coordinates of alpha carbon (for convenience in the KDtree class)
        float x() {return get_atom("CA")->x();}
        float y() {return get_atom("CA")->y();}
        float z() {return get_atom("CA")->z();}
        // distance to another residue CA atoms
        float distance (Residue&);
        // minimum distance between two residues
        float min_distance (Residue&);
        // bool operator to compare sizes (number of atoms)
        bool operator< (const Residue&) const;
        // function to count the number of hydrogen bonds between two residues
        // parameters are the other residue, the distance cutoff, and the angles
        // between donor-hydrogen-acceptor and donor-acceptor-antecedent
        vector<PROT::HydrogenBond> hbond(Residue *, float, float, float, bool);
        // function to find a salt bridge with another residue
        // the parameters are the other residue, the distance cutoff, and a verbose bool
        vector<SaltBridge> salt_bridge(Residue *, float, bool);
        // function to find hydrophobic interactions between two residues
        // the parameters are the other residue and a verbose bool
        vector<PROT::Hydrophobic> hydrophobic(Residue *, float, bool);
        // function to get the rotamers of this residue
        void set_rotamers();
        // function to get the rotamers of this residue (set them if m_rotamers is empty)
        vector<PROT::Residue> get_rotamers() {if (m_rotamers.size() == 0) {set_rotamers();} return m_rotamers;}
        size_t free_rotamers(vector<PROT::Residue>, vector<PROT::Residue*>);
        // function to get the neighbors of this residue
        vector<PROT::Residue*> get_intra_neighbors(PROT::Protein*, float);
        vector<PROT::Residue*> get_intra_neighbors(vector<PROT::Residue*>, float);
        vector<PROT::Residue*> get_inter_neighbors(vector<PROT::Protein*>, float);
        vector<PROT::Residue*> get_inter_neighbors_res(vector<PROT::Residue*>, float);
        // function to set the stability of the residue
        void set_stability(vector<PROT::Residue>, vector<PROT::Residue*>, bool, string&);
        // access the prestability of the residue
        bool prestable () {return m_prestable;}
        // access the bound stability of the residue
        bool bound_stable () {return m_bound_stable;}
        // function to check if two residues are clashing sterically
        bool clash (Residue*);
        bool heavy_clash (Residue*);
        bool heavy_side_chain_clash (Residue*);
        // remove side chains
        void remove_sidechain ();
        // calculate the rmsd between this residue and another
        float rmsd (PROT::Residue*);
    // End the Residue class definition
};

// Define the ResiduePtr class. This class holds a pointer to a Residue and
// provides access to all of it's public methods. 
class PROT::ResiduePtr {

    // The information stored in the class is private
    private:
        Residue * m_ptr;

    // A private methods to check that the pointer is non-null
    private:
        void check () const {
            if (m_ptr == 0) {
                string error = "Methods of the ResiduePtr class do not work "
                               "for Null pointers.\n";
                throw PANTZ_error (error);}}

    // The public interface of the class
    public:
        // Like the AtomPtr class, this class does NOT have a destructor,
        // because it is not the actual container of the information. It is just
        // a wrapper around a pointer to that data.
        // Class constructors
        ResiduePtr () {m_ptr = 0;}
        ResiduePtr (Residue * res) {m_ptr = res;}
        ResiduePtr (Residue& res) {m_ptr = &res;}
        ResiduePtr (const ResiduePtr& other) {m_ptr = other.m_ptr;}
        // Access to the pointer
        Residue * pointer () {check(); return m_ptr;}
        // Copy assignment
        void operator= (const ResiduePtr& other) {m_ptr = other.m_ptr;}
        void operator= (Residue * res) {m_ptr = res;}
        // Access to Residue methods
        size_t size () const {check(); return m_ptr->size();}
        string name () const {check(); return m_ptr->name();}
        char AA1 () const {check(); return m_ptr->AA1();}
        long number () const {check(); return m_ptr->number();}
        long internal_number () const {check(); return m_ptr->internal_number();}
        char insertion_code () const {check(); return m_ptr->insertion_code();}
        char protein () const {check(); return m_ptr->protein();}
        bool is_present () const {check(); return m_ptr->is_present();}
        size_t missing_atoms () const {check(); return m_ptr->missing_atoms();}
        coor phi () const {check(); return m_ptr->phi();}
        coor psi () const {check(); return m_ptr->psi();}
        coor omega () const {check(); return m_ptr->omega();}
        PROT::Atom * get_atom (const size_t i) {check(); return m_ptr->get_atom(i);}
        PROT::Atom * get_atom (const string& l) {check(); return m_ptr->get_atom(l);}
        PROT::AtomPtr operator[] (const size_t i) {
            check(); return m_ptr->operator[](i);}
        PROT::AtomPtr operator[] (const string& L) {
            check(); return m_ptr->operator[](L);}
        // Set number is defined with the set number method of the Residue class
        // in set.h to make sure that the default behaviors are the same
        void set_number (const long, const char, const bool);
        // Keep implementing other methods of the Residue class here
        void set_protein (const char L) {check(); m_ptr->set_protein(L);}
        long renumber_atoms (long n) {check(); return m_ptr->renumber_atoms(n);}
        long last_atom_number () const {check(); return m_ptr->last_atom_number();}
        // The string method has default behavior and is assigned in str.h
        string str (const bool);
        // The move functions also have default behaviors and are
        // defined with their Residue counterparts
        void move (const Matrix *, const char);
        void move (const Matrix *, const bool);
        void move (const Matrix&, const char);
        void move (const Matrix&, const bool);
        // Rotate functions can be implemented here
        void rotate (const Matrix * m) {check(); m_ptr->rotate(m);}
        void rotate (const Matrix& m) {check(); m_ptr->rotate(m);}
        // The select atoms methods also have default behaviors and are
        // implemented with their Residue counterparts
        void select_atoms (vector<Atom *>&, const string);
        void select_atoms (vector<AtomPtr>&, const string);
        // Center makes use of select atoms and has default behavior
        PROT::Matrix center (const string);
        // As does residue positioning
        void position (Matrix&, const bool);
        // The remaining methods can be implemented here
        double score () const {check(); return m_ptr->score();}
        bool is_amino_acid () const {check(); return m_ptr->is_amino_acid();}
        void for_charmm_histidine_fix () {
            check(); m_ptr->for_charmm_histidine_fix();}
        void from_charmm_histidine_fix() {
            check(); m_ptr->from_charmm_histidine_fix();}
        string rosetta_str (long& rn, long& an, const bool last) {
            check(); return m_ptr->rosetta_str(rn, an, last);}
        void update_atoms_after_Rosetta (const bool last) {
            check(); m_ptr->update_atoms_after_Rosetta(last);}
        Residue duplicate () const {check(); return m_ptr->duplicate();}

    // End the class definition
};

// Define a preprocessor variable to guarantee that the Residue's methods are
// being included from this file
#define Residue_Loading_Status 1

// (added by clay)
// include the types of interactions that can be found between residues
#include "HydrogenBond.h"
#include "SaltBridge.h"
#include "Hydrophobic.h"

// Include the methods of the Residue class
#include "Residue/initialize.h"
#include "Residue/clean_up.h"
#include "Residue/copy.h"
#include "Residue/set.h"
#include "Residue/check_atoms.h"
#include "Residue/move.h"
#include "Residue/rotate.h"
#include "Residue/load.h"
#include "Residue/constructors.h"
#include "Residue/dihedrals.h"
#include "Residue/AA1.h"
#include "Residue/operators.h"
#include "Residue/renumber_atoms.h"
#include "Residue/last_atom_number.h"
#include "Residue/str.h"
#include "Residue/select_atoms.h"
#include "Residue/center.h"
#include "Residue/position.h"
#include "Residue/score.h"
#include "Residue/histidine.h"
#include "Residue/rosetta.h"
#include "Residue/duplicate.h"
// added by clay
#include "Residue/distance.h"
#include "Residue/hbond.h"
#include "Residue/salt_bridge.h"
#include "Residue/hydrophobic.h"
#include "Residue/stability.h"
#include "Residue/neighbors.h"
#include "Residue/rotamers.h"
#include "Residue/clash.h"
#include "Residue/remove_sidechain.h"
#include "Residue/rename.h"
#include "Residue/rmsd.h"

// Include residue-based matrix allocation methods
#include "Matrix/allocate_residue.h"

// Undefine the loading status preprocessor variable
#undef Residue_Loading_Status

// End the header guard from the start of the file
#endif
