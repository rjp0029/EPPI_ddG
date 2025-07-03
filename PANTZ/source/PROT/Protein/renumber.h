/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * methods for renumbering the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Change the internal numbers of the Residues
long PROT::Protein::renumber_residues (long n) {
    // If there are Residues
    if (m_count > 0) {
        // Loop through them
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i].m_internal = n;
            ++n;}}
    return n;
}

// Change the numbers of the Atoms in the Protein
long PROT::Protein::renumber_atoms (long n) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            n = m_residues[i].renumber_atoms (n);}}
    return n;
}

