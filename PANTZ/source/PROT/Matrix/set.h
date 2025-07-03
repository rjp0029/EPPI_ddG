/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the set method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
void PROT::Matrix::set (const size_t r, const size_t c, const coor value) {
    // Error check the indices
    error_check(r, c);
    // Get a linear index
    size_t n = index(r, c);
    // Store the value
    m_values[n] = value;
}
