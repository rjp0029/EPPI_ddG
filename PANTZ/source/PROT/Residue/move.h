/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * move methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The private move function of the Residue class. It does not check that the
// provided matrix is valid.
void PROT::Residue::private_move (const Matrix * matrix, const char how) {
    // If there are Atoms, use their private move functions. The Residue class
    // is a friend of the Atom class, so this is permitted
    if (m_count > 0) {
        for (size_t i=0; i<m_count; ++i){m_atoms[i].private_move(matrix, how);}}
}

// Move the Residue, but error check the inputs, first
void PROT::Residue::move (const Matrix * matrix, const char how = '-') {
    // Check the matrix
    matrix->move_check();
    // Check the how character
    if ((how != '-') && (how != '+')) {
        string error = "'";
        error += how;
        error += "' is not a valid character for how to move Atoms.\n";
        throw PANTZ_error (error);}
    // Use the private move function
    private_move(matrix, how);
}

// Do the same things, but use bools instead of chars to indicate how to do
// things
void PROT::Residue::move (const Matrix * matrix, const bool subtract = true) {
    // Error check the matrix
    matrix->move_check();
    // Since it is a bool, it can only be true or false
    if (subtract) {private_move(matrix, '-');}
    else {private_move(matrix, '+');}
}

// Do the same things, but with matrices passed by reference
void PROT::Residue::move (const Matrix& matrix, const char how = '-') {
    move(&matrix, how);
}

void PROT::Residue::move (const Matrix& matrix, const bool subtract = true) {
    move(&matrix, subtract);
}

// The 4 equivalent methods in the ResiduePtr class
void PROT::ResiduePtr::move (const Matrix * matrix, const char how = '-') {
    check (); 
    m_ptr->move(matrix, how);
}

void PROT::ResiduePtr::move (const Matrix * matrix, const bool subtract = true) {
    check();
    m_ptr->move(matrix, subtract);
}

void PROT::ResiduePtr::move (const Matrix& matrix, const char how = '-') {
    check ();
    m_ptr->move(matrix, how);
}

void PROT::ResiduePtr::move (const Matrix& matrix, const bool subtract = true) {
    check();
    m_ptr->move(matrix, subtract);
}
