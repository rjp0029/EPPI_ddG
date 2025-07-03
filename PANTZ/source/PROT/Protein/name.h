/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the set_name method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Change the Protein's name
void PROT::Protein::set_name (const char L) {
    // Check that the name is valid
    CHECK::protein_name (L);
    // Modify the protein's attribute
    m_name = L;
    // If there are Residues, change them, too
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i].private_set_protein(L);}}
}
