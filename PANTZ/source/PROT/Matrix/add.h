/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the add method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
void PROT::Matrix::add (const size_t r, const size_t c, const coor value) {
    // Error check the indices
    error_check (r, c);
    // Calculate the linear index
    size_t n = index(r, c);
    // Add the value to the current value
    m_values[n] += value;
}
