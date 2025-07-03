/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the move methods of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Move the Atom without error checking the inputs
void PROT::Atom::private_move (const Matrix * matrix, const char how) {
    // The default behaviour is to subtract the coordinates from the matrix
    if (how == '-') {
        for(size_t i=0; i<AtomCoordinates; ++i) {
            m_coors[i] -= matrix->operator()(0, i);}}
    // If they aren't being subtracted, they should be added
    else {
        for(size_t i=0; i<AtomCoordinates; ++i) {
            m_coors[i] += matrix->operator()(0, i);}}
    return;
}

// Move the Atom using a pointer to a matrix and confirming that the Matrix is
// valid for the task
void PROT::Atom::move (const Matrix * matrix, const char how = '-') {
    matrix->move_check();
    // Error check that how is an acceptable value
    if ((how != '-') && (how != '+')) {
        string error = "'";
        error += how;
        error += "' is not an acceptable character for how to move an Atom.\n";
        throw PANTZ_error (error);}
    // Move the Atom
    private_move(matrix, how);
}

// Do the same, but use a Matrix by reference
void PROT::Atom::move (const Matrix& matrix, const char how = '-') {
    move (&matrix, how);
}

// Move functions for AtomPtrs
void PROT::AtomPtr::move (const Matrix * matrix, const char how = '-') {
    // Error check the pointer
    check();
    // Call the method of the atom
    m_ptr->move(matrix, how);
}

void PROT::AtomPtr::move (const Matrix& matrix, const char how = '-') {
    check ();
    m_ptr->move(matrix, how);
}
