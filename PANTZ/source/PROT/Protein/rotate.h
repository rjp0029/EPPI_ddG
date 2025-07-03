/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the rotation methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Rotate the protein without error checking the matrix
void PROT::Protein::private_rotate (const Matrix * mat) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_residues[i].private_rotate(mat);}}
}

// Rotate the protein with matrix error checking
void PROT::Protein::rotate (const Matrix * mat) {
    mat->rotate_check ();
    private_rotate(mat);
}

void PROT::Protein::rotate (const Matrix& matrix) {rotate(&matrix);}

