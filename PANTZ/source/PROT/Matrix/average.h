/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the code to calculate the 
 * average of a matrix either rowwise or columnwise. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// average the matrix rowwise or columnwise
PROT::Matrix PROT::Matrix::average (const bool rowwise = true) const {
    // if the matrix is empty, throw an error
    if (m_rows == 0 || m_columns == 0) {
        throw PANTZ_error ("It is not possible to calculate the average of an empty matrix.\n");}
    // make a new matrix
    Matrix output (1, rowwise ? m_columns : m_rows);
    // calculate the average
    for (size_t i=0; i<(rowwise ? m_columns : m_rows); ++i) {
        coor sum = 0;
        for (size_t j=0; j<(rowwise ? m_rows : m_columns); ++j) {
            sum += rowwise ? m_values[index(j, i)] : m_values[index(i, j)];
        }
        output.set(0, i, sum / (rowwise ? m_rows : m_columns));
    }
    return output;
}
