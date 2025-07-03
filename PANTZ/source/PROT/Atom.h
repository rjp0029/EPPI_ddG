/* Created by the Pantazes Lab at Auburn University.
 *
 * This file defines the Atom class for working with PDB-formatted information.
 * It is intended to be included by the Proteins.h header file and it includes a
 * number of other header files with the class methods. */

// Define a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_Atom_Guard
#define Proteins_Atom_Guard 1

// Confirm that the Proteins header file is being loaded
#ifndef Proteins_Loading_Status
#error Atom.h must be included from Proteins.h
#endif

// Confirm that the Check and Matrix files have been loaded
#ifndef Check_Guard
#error Atom.h must be inlcuded after Check.h
#endif

#ifndef Proteins_Matrix_Guard
#error Atom.h must be included after Matrix.h
#endif

// Define the Atom class
class PROT::Atom {

    // The Residue class is a friend of the Atom class so that it can directly
    // modify its parameters
    friend class PROT::Residue;
    // The PDB class also needs to be a friend
    friend class PROT::PDB;

    // The information stored in the class is private
    private:
        // The Atom's coordinates
        coor m_coors [AtomCoordinates];
        // ATOM or HETATM at the start of the line
        string m_type;
        // The atom's name
        string m_name;
        // The name of the atom's residue
        string m_residue;
        // The atom's chemical element
        string m_element;
        // The atom's charge (often missing, so stored as a string)
        string m_charge;
        // The atom's number
        long m_number;
        // The residue's number
        long m_residue_number;
        // The occupancy of the atom in the PDB file (x-ray crystallography
        // property)
        float m_occupancy;
        // The temperature factur
        float m_temperature;
        // The alternate location specifier for the atom's residue. This is
        // related to but not the same as the occupancy
        char m_alt;
        // The residue's insertion code - most relevant for antibody CDRs
        char m_insertion;
        // The atom's residue's protein's name
        char m_protein;
        // the atoms solvent accesible points (set by METHODS::set_sasa_points)
        vector<vector<float>> m_sasa_points;

    // Private methods to control class behavior
    private:
        // Assign default values to the class variables
        void initialize();
        // Copy information from another instance of the class
        void copy (const Atom *);
        void copy (const Atom& other) {copy(&other);}
        // Move the Atom without confirming the validity of the matrix for that
        // purpose
        void private_move (const Matrix *, const char);
        // Same idea for rotation
        void private_rotate (const Matrix *);

        // Energy calculations directly by this C++ code were in a previous
        // version, but have been removed because they're no longer directly
        // used. However, while the individual energy calculations are simple,
        // determining when they should be ignored is not. This function is for
        // that purpose. It is commented out here and it's header file will not
        // be included in the code, but they're being left in place in case
        // those features are ever re-added to the code
        // int exclusions (const Atom *) const;

    // The public interface of the class
    public:
        // The default class constructor
        Atom () {initialize();}
        // The standard constructor
        Atom (const string&);
        // Copy construction and assignment
        Atom (const Atom& other) {copy(other);}
        void operator= (const Atom& other) {copy(other);}
        // Access to class attributes
        size_t size () const {return AtomCoordinates;}
        string type () const {return m_type;}
        string name () const {return m_name;}
        string residue () const {return m_residue;}
        string element () const {return m_element;}
        string charge () const {return m_charge;}
        long number () const {return m_number;}
        long residue_number () const {return m_residue_number;}
        float occupancy () const {return m_occupancy;}
        float temperature () const {return m_temperature;}
        char alternative_location () const {return m_alt;}
        char insertion_code () const {return m_insertion;}
        char protein () const {return m_protein;}
        // Access to the atom's coordinates
        coor operator[] (const size_t) const;
        coor operator[] (const char) const;
        // Create a PDB formatted string of the Atom's data
        string str () const;
        // Create a PDB formatted string of the Atom's data that is appropriate
        // for use in Rosetta
        string rosetta_str () const;
        // Whether or not the Atom is a backbone atom
        bool is_backbone_atom () const {return CHECK::is_backbone_atom(m_name);}
        // Whether or not the atom is a hydrogen
        bool is_hydrogen () const;
        // The chemical element of the atom (the first non-digit character of
        // it's name). This function is implemented in the is_hydrogen header
        // file
        char determine_element () const;
        // Update an Atom's name for use in Rosetta
        void update_name_for_Rosetta (const bool);
        // And after Rosetta
        void update_name_after_Rosetta (const bool);
        // Move the Atom, checking the matrix for errors first
        void move (const Matrix *, const char);
        void move (const Matrix&, const char);
        void rotate (const Matrix *);
        void rotate (const Matrix& matrix) {rotate(&matrix);}
        // Calculate the distance between two atoms
        coor calculate_distance (const Atom *, const bool) const;
        coor calculate_distance (const Atom&, const bool) const;
        // change the residue number of the atom
        void renumber_residue (const long n) {m_residue_number = n;}
        // x y and z coordinates (for convenience in the KDtree class)
        float x() const {return m_coors[0];} 
        float y() const {return m_coors[1];} 
        float z() const {return m_coors[2];}
        // get the distance to another atom
        float distance (const Atom& other) const {return calculate_distance(other, false);}
        // the sigma value for the atom
        float lj_sigma() const;
        // the sphere around the atom (for interaction calculations)
        vector<vector<float>> fibonacci_sphere (size_t N, float probe_radius);
        // set sasa points to the atom
        void add_sasa_point (vector<float> point) {m_sasa_points.push_back(point);}
        // check if atoms are sterically clashing
        bool clash (Atom*, float probe_radius = 0);

    // End the class definition
};

// The AtomPtr class is a wrapper class that just holds a pointer to an Atom. It
// has the same public list of methods as an Atom (and related ones with other
// AtomPtrs) to allow for simpler coding elsewhere in PANTZ
class PROT::AtomPtr {

    // The information stored in the class is private
    private:
        Atom * m_ptr;

    // There is a single private method that confirms the pointer is non-Null
    private:
        void check () const {
            if (m_ptr == 0) {
                string error = "Methods of the AtomPtr class do not work for "
                               "Null pointers.\n";
                throw PANTZ_error (error);}}

    // The public interface of the class
    public:
        // The class does NOT have a destructor. It is not the master copy of an
        // Atom, just a version that can be passed by reference to it. 
        // Class constructors
        AtomPtr () {m_ptr = 0;}
        AtomPtr (Atom * atom) {m_ptr = atom;}
        AtomPtr (Atom& atom) {m_ptr = &atom;}
        AtomPtr (const AtomPtr& other) {m_ptr = other.m_ptr;}
        // Access to the pointer of the class
        Atom * pointer () {check(); return m_ptr;}
        // Copy assignment to an AtomPtr
        void operator= (AtomPtr& other) {m_ptr = other.m_ptr;}
        // Copy assignment to Atoms
        void operator= (Atom * atom) {m_ptr = atom;}
        // Access to the Atom's functions
        size_t size () const {check(); return m_ptr->size();}
        string type () const {check(); return m_ptr->type();}
        string name () const {check(); return m_ptr->name();}
        string residue () const {check(); return m_ptr->residue();}
        string element () const {check(); return m_ptr->element();}
        string charge () const {check(); return m_ptr->charge();}
        long number() const {check(); return m_ptr->number();}
        long residue_number () const {check(); return m_ptr->residue_number();}
        float occupancy () const {check(); return m_ptr->occupancy();}
        float temperature () const {check(); return m_ptr->temperature();}
        char alternative_location () const {check(); return m_ptr->alternative_location();}
        char insertion_code () const {check(); return m_ptr->insertion_code();}
        char protein () const {check(); return m_ptr->protein();}
        coor operator[] (const size_t i) const {check(); return m_ptr->operator[](i);}
        coor operator[] (const char l) const {check(); return m_ptr->operator[](l);}
        string str () const {check(); return m_ptr->str();}
        string rosetta_str () const {check(); return m_ptr->rosetta_str();}
        bool is_backbone_atom () const {check(); return m_ptr->is_backbone_atom();}
        bool is_hydrogen () const {check(); return m_ptr->is_hydrogen();}
        char determine_element () const {check(); return m_ptr->determine_element();}
        void update_name_for_Rosetta (const bool last = false) const {
            check (); m_ptr->update_name_for_Rosetta(last);}
        void update_name_after_Rosetta (const bool last = false) const {
            check (); m_ptr->update_name_after_Rosetta(last);}
        // The move, rotate, and calculate distance functions are defined in the
        // corresponding atom header files to make sure that the default
        // behaviors remain consistent.
        void move (const Matrix *, const char);
        void move (const Matrix&, const char);
        void rotate (const Matrix *);
        void rotate (const Matrix&);
        coor calculate_distance (const Atom *, const bool) const;
        coor calculate_distance (const Atom&, const bool) const;
        coor calculate_distance (const AtomPtr&, const bool) const;

    // End the class definition
};

// Define a pre-processor variable to guarantee that Atom methods are loaded
// from this file
#define Atom_Loading_Status

// Include the Atom's header files
#include "Atom/initialize.h"
#include "Atom/copy.h"
#include "Atom/move.h"
#include "Atom/rotate.h"
#include "Atom/constructor.h"
#include "Atom/operator.h"
#include "Atom/str.h"
#include "Atom/is_hydrogen.h"
// added by clay
#include "Atom/calculate_distance.h"
#include "Atom/update_name.h"
#include "Atom/fibonacci_sphere.h"
#include "Atom/clash.h"

// Include the atom allocation methods of the Matrix class
#include "Matrix/allocate_atom.h"

// Include the calculate dihedral functions
#include "Atom/calculate_dihedral.h"

// Undefine the atom loading status
#undef Atom_Loading_Status

// End the header guard from near the start of the file
#endif
