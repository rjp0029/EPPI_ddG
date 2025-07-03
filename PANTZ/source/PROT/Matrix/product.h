/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the dot and cross product methods of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Dot product
PROT::Matrix PROT::Matrix::dotProduct (const Matrix * other) const {
    // Validate that the row and column information matches up
    if (m_columns != other->m_rows) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        stringstream c3; c3 << other->m_rows;
        stringstream c4; c4 << other->m_columns;
        string error = "It is not possible to calculate the dot product of a "
                     + c1.str() + "x" + c2.str() + " matrix and a "
                     + c3.str() + "x" + c4.str() + " matrix.\n";
        throw PANTZ_error (error);}
    // Make the new matrix
    Matrix output (m_rows, other->m_columns);
    // A counter
    size_t n = 0;
    // Loop through the rows of this matrix
    for(size_t i=0; i<m_rows; ++i) {
        // Loop through the columns of the other matrix
        for(size_t j=0; j<other->m_columns; ++j) {
            // The number of terms in the sum
            for(size_t k=0; k<m_columns; ++k) {
                size_t N1 = index(i, k);
                size_t N2 = other->index(k, j);
                output.m_values[n] += m_values[N1] * other->m_values[N2];}
            ++n;}}
    return output;
}

// Dot product using a matrix passed by reference
PROT::Matrix PROT::Matrix::dotProduct (const Matrix& other) const {
    return dotProduct(&other);
}

// Cross product
PROT::Matrix PROT::Matrix::crossProduct (const Matrix * O) const {
    // Current calculations only support 1x3 calculations. Check the dimensions
    // of both
    bool flag1 = false;
    if (((m_rows == 1) && (m_columns == AtomCoordinates)) ||
        ((m_rows == AtomCoordinates) && (m_columns == 1))) {flag1 = true;}
    bool flag2 = false;
    if (((O->m_rows == 1) && (O->m_columns == AtomCoordinates)) ||
        ((O->m_rows == AtomCoordinates) && (O->m_columns == 1))) {flag2 = true;}
    // If either flag is false, the dimensions don't match
    if ((!flag1) || (!flag2)) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Currently, cross products can only be calculated "
                       "between 1x" + c1.str() + " vectors.\n";
        throw PANTZ_error (error);}
    // Create an output matrix
    Matrix output (1, AtomCoordinates);
    // Calculate the values
    output.m_values[0] = m_values[1]*O->m_values[2] - m_values[2]*O->m_values[1];
    output.m_values[1] = m_values[2]*O->m_values[0] - m_values[0]*O->m_values[2];
    output.m_values[2] = m_values[0]*O->m_values[1] - m_values[1]*O->m_values[0];
    return output;
}

PROT::Matrix PROT::Matrix::crossProduct (const Matrix& other) const {
    return crossProduct(&other);
}
