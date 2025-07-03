/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to remove side chain atoms. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// remove the side chain atoms
void PROT::Residue::remove_sidechain () {
    size_t bbcount = 0;
    // iterate through the atoms
    for (size_t i = 0; i < m_count; i++) {
        // if the atom is not a backbone atom
        if (m_atoms[i].is_backbone_atom()) {
            bbcount++;
        }
    }
    // create the new array of atoms
    PROT::Atom * m_atoms_new = new PROT::Atom [bbcount];
    // iterate through the atoms
    size_t j = 0;
    for (size_t i = 0; i < m_count; i++) {
        // if the atom is not a backbone atom
        if (m_atoms[i].is_backbone_atom()) {
            // copy the atom to the new array
            m_atoms_new[j] = m_atoms[i];
            j++;
        }
    }
    // delete the old array
    delete [] m_atoms;
    // set the new array
    m_atoms = m_atoms_new;
    // set the new count
    m_count = bbcount;
    return;
}