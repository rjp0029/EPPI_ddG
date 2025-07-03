/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * renaming method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Renumber the Atoms in the Residue
void PROT::Residue::rename(const string& new_name) {
    m_name = new_name;
    // iterate through the atoms and set m_residue to new name
    for (size_t i = 0; i < m_count; i++) {
        m_atoms[i].m_residue = new_name;
    }
    // return the number of atoms
}
