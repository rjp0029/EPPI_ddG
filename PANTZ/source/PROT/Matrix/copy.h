/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the copy method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the pointer-based method
void PROT::Matrix::copy (const Matrix * other) {
    // Delete any existing information
    clean_up();
    // Store the dimensions
    m_rows = other->m_rows;
    m_columns = other->m_columns;
    // Calculate the total number of entries
    size_t n = m_rows * m_columns;
    // If there are values
    if (n > 0) {
        // Allocate memory
        m_values = new coor [n];
        // Copy them
        for(size_t i=0; i<n; ++i) {
            m_values[i] = other->m_values[i];}}
}

