/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the duplicate method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Create a duplicated copy of a Protein
PROT::Protein PROT::Protein::duplicate () const {
    // Create a new protein
    Protein output;
    // Copy the name
    output.m_name = m_name;
    // Copy the number of residues
    output.m_count = m_count;
    // If there are residues
    if (m_count > 0) {
        // Allocate memory
        output.m_residues = new Residue [m_count];
        // Copy the residues
        for(size_t i=0; i<m_count; ++i) {
            output.m_residues[i] = m_residues[i];}}
    return output;
}
