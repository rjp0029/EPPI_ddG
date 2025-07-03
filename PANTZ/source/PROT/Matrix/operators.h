/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the "operator" methods of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Access to values by indexes
PROT::coor PROT::Matrix::operator() (const size_t r, const size_t c) const {
    // Error check the indices provided
    error_check (r, c);
    // Convert them to a linear index
    size_t n = index(r, c);
    // Return the value
    return m_values[n];
}

// Matrix subtraction
PROT::Matrix PROT::Matrix::operator- (const Matrix& other) const {
    // Make sure the dimensions match
    if ((m_rows != other.m_rows) || (m_columns != other.m_columns)) {
        string error = "Matrix dimensions must match for subtraction to be a "
                       "valid operation.\n";
        throw PANTZ_error (error);}
    // Make sure there are values
    if ((m_rows == 0) || (m_columns == 0)) {
        string error = "Subtraction is not a valid operation for an empty "
                       "matrix.\n";
        throw PANTZ_error (error);}
    // Make a new matrix
    Matrix output (m_rows, m_columns);
    // Assign the values
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        output.m_values[i] = m_values[i] - other.m_values[i];}
    return output;
}

// Matrix addition
PROT::Matrix PROT::Matrix::operator+ (const Matrix& other) const {
    if ((m_rows != other.m_rows) || (m_columns != other.m_columns)) {
        string error = "Matrix dimensions must match for addition to be a "
                       "valid operation.\n";
        throw PANTZ_error (error);}
    if ((m_rows == 0) || (m_columns == 0)) {
        string error = "Addition is not a valid operation for an empty "
                       "Matrix.\n";
        throw PANTZ_error (error);}
    Matrix output (m_rows, m_columns);
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        output.m_values[i] = m_values[i] + other.m_values[i];}
    return output;
}
