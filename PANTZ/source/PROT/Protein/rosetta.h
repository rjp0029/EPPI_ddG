/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the Rosetta-related methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Make a string of information formatted for use in Rosetta
string PROT::Protein::rosetta_str (long& resNum, long& atomNum) {
    // Store the output here
    string output = "";
    // Call the methods of the Residues
    for(size_t i=0; i<m_count; ++i) {
        output.append(m_residues[i].rosetta_str(resNum, atomNum, i==m_count-1));}
    // Return the string
    return output;
}

// Update atom naming information after Rosetta calculations
void PROT::Protein::update_atoms_after_Rosetta () {
    for(size_t i=0; i<m_count; ++i) {
        m_residues[i].update_atoms_after_Rosetta(i==m_count-1);}
}
