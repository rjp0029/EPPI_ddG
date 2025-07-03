/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the initialize method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
void PROT::Matrix::initialize () {
    // Set the numbers of rows and columns to 0
    m_rows = 0;
    m_columns = 0;
    // Set the values array to empty
    m_values = 0;
}
