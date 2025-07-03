/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the rotate methods of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Rotate the Atom's position without error checking the matrix
void PROT::Atom::private_rotate (const Matrix * matrix) {
    // Create an array of the new coordinates
    coor newCoors [AtomCoordinates];
    // Calculate the new coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Set the new coordinate value to 0
        newCoors[i] = 0;
        // Loop through the current coordinates
        for(size_t j=0; j<AtomCoordinates; ++j) {
            // Add the appropriate product value to the new coordinate
            newCoors[i] += (matrix->operator()(i, j) * m_coors[j]);}}
    // Store the new coordinates as the Atom's coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = newCoors[i];}
}

// Rotate an Atom, but only after checking that the provide matrix is
// appropriate for that calculation.
void PROT::Atom::rotate (const Matrix * matrix) {
    matrix->rotate_check ();
    private_rotate(matrix);
}

// Rotation functions for AtomPtrs
void PROT::AtomPtr::rotate (const Matrix * matrix) {
    // Error check the atom pointer
    check ();
    // Call the rotate method of the atom object
    m_ptr->rotate(matrix);
}

void PROT::AtomPtr::rotate (const Matrix& matrix) {
    check ();
    m_ptr->rotate(matrix);
}
