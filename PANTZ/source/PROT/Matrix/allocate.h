/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the allocate methods of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the methods

// Allocate space for the matrix to be r x c
void PROT::Matrix::allocate (const size_t r, const size_t c) {
    // Error check the number of rows and columns
    if (r < 1) {
        stringstream c1; c1 << r;
        string error = "A Matrix must have a positive number of rows.\n"
                     + c1.str() + " is not acceptable.\n";
        throw PANTZ_error (error);}
    else if (c < 1) {
        stringstream c1; c1 << c;
        string error = "A Matrix must have a positive number of columns.\n"
                     + c1.str() + " is not acceptable.\n";
        throw PANTZ_error (error);}
    // Delete current information
    clean_up();
    // Assign the values and set each entry to 0
    m_rows = r;
    m_columns = c;
    size_t n = m_rows * m_columns;
    m_values = new coor [n];
    for(size_t i=0; i<n; ++i) {m_values[i] = 0;}
}

// Set up a Matrix for a rotation calculation using Rodriguez's rotation
// formula
void PROT::Matrix::allocate (const coor angle, const coor vector []) {
    // Set up the Matrix to have appropriate dimensions
    clean_up();
    m_rows = AtomCoordinates;
    m_columns = AtomCoordinates;
    m_values = new coor [m_rows * m_columns];
    // Calculate some trig values that are needed
    coor c = cos(angle);
    coor s = sin(angle);
    coor v = 1 - c;
    // Store values in the matrix
    m_values[0] = c + v * vector[0] * vector[0];
    m_values[1] = -s*vector[2] + v * vector[0] * vector[1];
    m_values[2] = s*vector[1] + v * vector[0] * vector[2];
    m_values[3] = s*vector[2] + v * vector[1] * vector[0];
    m_values[4] = c + v * vector[1] * vector[1];
    m_values[5] = -s*vector[0] + v * vector[1] * vector[2];
    m_values[6] = -s*vector[1] + v * vector[2] * vector[0];
    m_values[7] = s*vector[0] + v * vector[2] * vector[1];
    m_values[8] = c + v * vector[2] * vector[2];
    return;
}

