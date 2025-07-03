/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the deallocation method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Delete dynamically allocated memory
void PROT::Protein::clean_up () {
    if (m_residues != 0) {delete[] m_residues; m_residues = 0;}
}
