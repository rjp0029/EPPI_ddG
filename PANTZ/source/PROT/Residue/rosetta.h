/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * Rosetta-specific methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Create a string of formatted text for the Residue's non-hydrogen atoms for
// use in Rosetta
string PROT::Residue::rosetta_str (long& resNum, long& atomNum, 
                                   const bool lastRes) {
    // Store the output string here
    string output = "";
    // Make sure histidine is HIS
    from_charmm_histidine_fix ();
    // Loop through the Residue's atoms
    for(size_t i=0; i<m_count; ++i) {
        // // If the atom is a hydrogen, skip it
        // if (m_atoms[i].is_hydrogen()) {continue;}
        // Update it's name for Rosetta
        m_atoms[i].update_name_for_Rosetta(lastRes);
        // Update the atom's number
        m_atoms[i].m_number = atomNum;
        atomNum++;
        // Use the proper residue numbering information
        if (resNum > 0) {
            m_atoms[i].m_residue_number = resNum;
            m_atoms[i].m_insertion = ' ';}
        else {
            m_atoms[i].m_residue_number = m_number;
            m_atoms[i].m_insertion = m_insertion;}
        // Add it's string to the output string
        output.append(m_atoms[i].rosetta_str());}
    // If appropriate, increment the residue number
    if (resNum > 0) {resNum++;}
    // Return the string
    return output;
}

// Update atom names after Rosetta
void PROT::Residue::update_atoms_after_Rosetta (const bool lastRes) {
    // Use the methods of the atoms
    for(size_t i=0; i<m_count; ++i) {
        m_atoms[i].update_name_after_Rosetta(lastRes);}
}
