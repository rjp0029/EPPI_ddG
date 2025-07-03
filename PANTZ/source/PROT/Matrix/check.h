/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the checking methods of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the indexing error check
void PROT::Matrix::error_check (const size_t r, const size_t c) const {
    // Check the provided row index. Note that size_t values are unsigned and
    // therefore the possibility of a value being negative does not need to be
    // checked
    if (r >= m_rows) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << r;
        string error = c2.str() + " is not a valid row index in a Matrix with "
                     + c1.str() + " rows.\n";
        throw PANTZ_error (error);}
    // Do the same for the column index
    if (c >= m_columns) {
        stringstream c1; c1 << m_columns;
        stringstream c2; c2 << c;
        string error = c2.str() + " is not a valid column index in a Matrix with "
                     + c1.str() + " columns.\n";
        throw PANTZ_error (error);}
}

// Whether the matrix is valid for moving an Atom
void PROT::Matrix::move_check () const {
    // Validate that it is a 1x3 matrix
    if (((m_rows != 1) || (m_columns != AtomCoordinates)) &&
        ((m_rows != AtomCoordinates) || (m_columns != 1))) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Only a 1x" + c1.str() + " vector can be used in "
                       "movement calculations.\n";
        throw PANTZ_error (error);}
}

// Whether the matrix is a valid rotation matrix
void PROT::Matrix::rotate_check () const {
    if ((m_rows != AtomCoordinates) || (m_columns != AtomCoordinates)) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Only a " + c1.str() + "x" + c1.str() + " Matrix can be "
                       "used in rotation calculations.\n";
        throw PANTZ_error (error);}
}
