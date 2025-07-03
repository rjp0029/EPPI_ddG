/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the score method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Calculate a completeness score for the Protein. This is only useful in the
// context of PDB files
double PROT::Protein::score () const {
    // Store the score here
    double value = 0.0;
    // If there are no residues, be done
    if (m_count == 0) {return value;}
    // Loop through the Residues
    for(size_t i=0; i<m_count; ++i) {
        // Add the score for the residue to the score for the protein
        value += m_residues[i].score();}
    // Return the value divided by the number of residues
    return value / ((double) m_count);
}
