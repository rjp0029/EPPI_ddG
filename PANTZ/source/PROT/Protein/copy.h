/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the copy method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Copy the information from another Protein into this one
void PROT::Protein::copy (const Protein * other) {
    // Delete any existing information
    clean_up();
    // Copy the name
    m_name = other->m_name;
    // Copy the number of Residues
    m_count = other->m_count;
    // If there are Residues, allocate memory and copy them
    if (m_count > 0) {
        m_residues = new Residue [m_count];
        for(size_t i=0; i<m_count; ++i) {m_residues[i] = other->m_residues[i];}}
}
