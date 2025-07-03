/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the string method of the class. */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

// Implement the method
string PROT::Matrix::str () const {
    // Store the output here
    string output = "";
    // If there are no values
    if (m_values == 0) {return output;}
    // Loop through the positions
    for(size_t i=0; i<m_rows; ++i) {
        for(size_t j=0; j<m_columns; ++j) {
            // Make a fixed precision representation
            stringstream c1; 
            c1 << fixed << setprecision(3) << m_values[index(i, j)];
            Text::rjust_insert(output, c1.str(), 12, ' ');}
        output += "\n";}
    return output;
}
