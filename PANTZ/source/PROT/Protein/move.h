/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the move methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// A private method that moves the protein without error checking either the
// matrix or the how character
void PROT::Protein::private_move (const Matrix * mat, const char how) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_residues[i].private_move(mat, how);}}
}

// Using a Matrix, move the Protein through space. Include error checking that
// the matrix is appropriate for this task
void PROT::Protein::move (const Matrix * mat, const char how = '-') {
    // Error check the matrix
    mat->move_check();
    // use the private move function
    private_move (mat, how);
}

// Use a boolean value instead of a character
void PROT::Protein::move (const Matrix * mat, const bool subtract = true) {
    if (subtract) {move(mat, '-');}
    else {move(mat, '+');}
}

// Do the same things with reference matrices instead of matrix pointers
void PROT::Protein::move (const Matrix& matrix, const char how = '-') {
    move(&matrix, how);
}

void PROT::Protein::move (const Matrix& matrix, const bool subtract = true) {
    move(&matrix, subtract);
}

