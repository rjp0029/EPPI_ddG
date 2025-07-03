/* Created by clay at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the formula to get the angle between two 1,3 matrices. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// get the angle between two 1,3 matrices
float PROT::Matrix::angle (const Matrix * other) const {
    // Validate that the row and column information matches up
    if ((m_rows != 1) || (m_columns != AtomCoordinates) ||
        (other->m_rows != 1) || (other->m_columns != AtomCoordinates)) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        stringstream c3; c3 << other->m_rows;
        stringstream c4; c4 << other->m_columns;
        string error = "It is not possible to calculate the angle between a "
                     + c1.str() + "x" + c2.str() + " matrix and a "
                     + c3.str() + "x" + c4.str() + " matrix.\n";
        throw PANTZ_error (error);}
    // Calculate the dot product of the two matrices
    float dot = 0;
    for (size_t i=0; i<AtomCoordinates; ++i) {
        dot += m_values[i] * other->m_values[i];}
    // Calculate the magnitudes of the two matrices
    float mag1 = 0;
    float mag2 = 0;
    for (size_t i=0; i<AtomCoordinates; ++i) {
        mag1 += m_values[i] * m_values[i];
        mag2 += other->m_values[i] * other->m_values[i];}
    mag1 = sqrt(mag1);
    mag2 = sqrt(mag2);
    // Calculate the angle
    return acos(dot / (mag1 * mag2));
}

    