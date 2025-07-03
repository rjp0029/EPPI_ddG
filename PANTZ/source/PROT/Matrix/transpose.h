/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the transpose method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
PROT::Matrix PROT::Matrix::transpose () const {
    // Validate that the matrix has information
    if ((m_rows == 0) || (m_columns == 0)) {
        string error = "An empty matrix cannot be transposed.\n";
        throw PANTZ_error (error);}
    // Make the matrix
    Matrix output (m_columns, m_rows);
    // Assign the values
    for(size_t i=0; i<m_rows; ++i) {
        for (size_t j=0; j<m_columns; ++j) {
            output.m_values[output.index(j, i)] = m_values[index(i, j)];}}
    return output;
}
