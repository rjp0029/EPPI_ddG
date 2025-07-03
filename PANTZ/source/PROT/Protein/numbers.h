/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * Protein class methods related to atom and residue numbers */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Calculate the number of Atoms in a Protein
size_t PROT::Protein::number_of_atoms () const {
    // Store the answer here
    size_t answer = 0;
    // If there are Residues
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {answer += m_residues[i].size();}}
    return answer;
}

// The number of the last Atom in the Protein
long PROT::Protein::last_atom_number () const {
    // If there are residues
    if (m_count > 0) {
        // Return the value from the last residue
        return m_residues[m_count-1].last_atom_number();}
    // Otherwise return 0
    return 0;
}

// The INTERNAL number of the last Residue in the Protein
long PROT::Protein::last_residue_number () const {
    if (m_count > 0) {return m_residues[m_count-1].internal_number();}
    return 0;
}
