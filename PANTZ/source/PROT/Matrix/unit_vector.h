/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the make unit vector method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
void PROT::Matrix::make_unit_vector () {
    // Validate that the matrix is a vector
    if ((m_rows != 1) && (m_columns != 1)) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        string error = "A " + c1.str() + "x" + c2.str() + " matrix is not a "
                       "vector.\n";
        throw PANTZ_error (error);}
    // Calculate the magnitude of the vector
    coor magnitude = 0.0;
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        magnitude += pow(m_values[i], 2);}
    magnitude = sqrt(magnitude);
    // Divide each value by the magnitude
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_values[i] /= magnitude;}
}
