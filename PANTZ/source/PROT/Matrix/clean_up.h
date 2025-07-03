/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the clean_up method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
void PROT::Matrix::clean_up () {
    // If there is dynamically allocated memory
    if (m_values != 0) {
        // Delete it
        delete[] m_values;
        // Set the pointer to 0
        m_values = 0;}
}
