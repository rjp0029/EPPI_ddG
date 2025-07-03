/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the declaration of the Matrix class for PDB calculations.
 * It also includes the header files where the methods of the class are
 * implemented. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_Matrix_Guard
#define Proteins_Matrix_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error Matrix.h must be included by Proteins.h
#endif

// Define the Matrix class
class PROT::Matrix {

    // The information stored in the class is private
    private:
        // The numbers of rows and columns
        size_t m_rows, m_columns;
        // The values in the matrix
        coor * m_values;

    // Private methods that control class behavior
    private:
        // Assign default values to class attributes
        void initialize ();
        // Deallocate dynamically assigned memory
        void clean_up ();
        // Copy information from another class instance
        void copy (const Matrix *);
        void copy (const Matrix& other) {copy(&other);}
        // Error check that row and column indices are valid
        void error_check (const size_t, const size_t) const;
        // Convert row and column information into a linear index
        size_t index (const size_t r, const size_t c) const {
            return c + (r*m_columns);}

    // The public interface of the class
    public:
        // The class destructor
        ~Matrix () {clean_up();}
        // The default class constructor
        Matrix () {initialize();}
        // Different methods of assigning values to a Matrix
        // By dimension
        void allocate (const size_t, const size_t);
        // Using Rodriguez's Rotation Formula
        void allocate (const coor, const coor []);
        // From an Atom
        void allocate (const Atom *);
        void allocate (const Atom& atom) {allocate(&atom);}
        void allocate (AtomPtr&);
        // Allocate a Matrix from a Residue
        void allocate (Residue *, const string);
        void allocate (Residue&, const string);
        void allocate (ResiduePtr&, const string);
        // From a vector of Atoms
        void allocate (const vector<Atom *>&);
        void allocate (vector<AtomPtr>&);
        // Constructors that use those allocation methods
        Matrix (const size_t r, const size_t c) {initialize(); allocate(r, c);}
        Matrix (const coor a, const coor v []) {initialize(); allocate(a, v);}
        Matrix (const Atom * atom) {initialize(); allocate(atom);}
        Matrix (const Atom& atom) {initialize(); allocate(atom);}
        Matrix (AtomPtr& atom) {initialize(); allocate(atom);}
        Matrix (const vector<Atom *>& atoms) {initialize(); allocate(atoms);}
        Matrix (vector<AtomPtr>& atoms) {initialize(); allocate(atoms);}
        Matrix (Residue *, const string);
        Matrix (Residue&, const string);
        Matrix (ResiduePtr&, const string);
        // Copy construction and assignment of a Matrix
        Matrix (const Matrix * other) {initialize(); copy(other);}
        Matrix (const Matrix& other) {initialize(); copy(other);}
        void operator= (const Matrix& other) {copy(other);}
        // Access to information in the matrix
        size_t rows () const {return m_rows;}
        size_t columns () const {return m_columns;}
        coor operator () (const size_t, const size_t) const;
        // Set values in a matrix
        void set (const size_t, const size_t, const coor);
        // Addition and subtraction of two matrices
        Matrix operator- (const Matrix&) const;
        Matrix operator+ (const Matrix&) const;
        // Addition to a single element
        void add (const size_t, const size_t, const coor);
        // Multiply two matrices together
        Matrix dotProduct (const Matrix *) const;
        Matrix dotProduct (const Matrix&) const;
        Matrix crossProduct (const Matrix *) const;
        Matrix crossProduct (const Matrix&) const;
        // Convert a vector into a unit vector
        void make_unit_vector ();
        // Transpose a matrix
        Matrix transpose () const;
        // Confirm that a matrix can be used to move atoms
        void move_check () const;
        // Confirm that a matrix can be used to rotate atoms
        void rotate_check () const;
        // A string representation of the Matrix
        string str () const;
        // Calculate the rotation matrix to rotate this matrix to the reference
        Matrix rotation_matrix (const Matrix*) const;
        // the average of the values in the matrix (can be rowwise or columnwise)
        Matrix average (const bool) const;
        // get the angle between two (1,3) matrices
        float angle (const Matrix*) const;
    // End the class definition
};

// Define a preprocessor variable to guarantee that the Matrix methods are
// included here and only here
#define Matrix_Loading_Status 1

// Include the files that implement class methods
#include "Matrix/initialize.h"
#include "Matrix/clean_up.h"
#include "Matrix/copy.h"
#include "Matrix/check.h"
#include "Matrix/allocate.h"
#include "Matrix/operators.h"
#include "Matrix/set.h"
#include "Matrix/add.h"
#include "Matrix/product.h"
#include "Matrix/unit_vector.h"
#include "Matrix/transpose.h"
#include "Matrix/str.h"
#include "Matrix/rotation_matrix.h"
#include "Matrix/average.h"
#include "Matrix/angle.h"

// Undefine the matrix loading status variable
#undef Matrix_Loading_Status

// End the header guard from the start of the file
#endif
