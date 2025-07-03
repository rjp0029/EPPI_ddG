/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * rotate methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The private rotate method of the class does not error check the provided
// matrix
void PROT::Residue::private_rotate (const Matrix * matrix) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            m_atoms[i].private_rotate(matrix);}}
}

// The public method does error check the matrix
void PROT::Residue::rotate (const Matrix * matrix) {
    matrix->rotate_check();
    private_rotate (matrix);
}

// Another public method using a reference matrix instead of a pointer
void PROT::Residue::rotate (const Matrix& matrix) {
    rotate(&matrix);
}
