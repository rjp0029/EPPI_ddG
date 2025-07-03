/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * methods to determine steric clashes between residues */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// determine if the atoms of this residue clash with the atoms of another residue
bool PROT::Residue::clash(PROT::Residue* other) {
    // iterate through m_atoms backwards and check for clashes
    for (size_t i = m_count; i > 0; i--) {
        // get the atom
        PROT::Atom * atom = &m_atoms[i-1];
        // iterate through the other residue's atoms
        for (size_t j = other->m_count; j > 0; j--) {
            // get the other atom
            PROT::Atom * other_atom = &other->m_atoms[j-1];
            // check for a clash
            if (atom->clash(other_atom)) {
                return true;
            }
        }
    }
    return false;
}

// determine if the atoms of this residue clash with the atoms of another residue
bool PROT::Residue::heavy_clash(PROT::Residue* other) {
    // iterate through m_atoms backwards and check for clashes
    for (size_t i = m_count; i > 0; i--) {
        // get the atom
        PROT::Atom * atom = &m_atoms[i-1];
        // if the atom is a hydrogen, skip it
        if (atom->m_element == "H") {
            continue;
        }
        // iterate through the other residue's atoms
        for (size_t j = other->m_count; j > 0; j--) {
            // get the other atom
            PROT::Atom * other_atom = &other->m_atoms[j-1];
            // if the other atom is a hydrogen, skip it
            if (other_atom->m_element == "H") {
                continue;
            }
            // check for a clash
            if (atom->clash(other_atom)) {
                return true;
            }
        }
    }
    return false;
}

// determine if the side chain atoms of this residue clash with the atoms of another residue
bool PROT::Residue::heavy_side_chain_clash(PROT::Residue* other) {
    // iterate through m_atoms backwards and check for clashes
    for (size_t i = m_count; i > 0; i--) {
        // get the atom
        PROT::Atom * atom = &m_atoms[i-1];
        
        // if the atom is a hydrogen, skip it
        if (atom->m_element == "H" or atom->is_backbone_atom()) {
            continue;
        }
        // iterate through the other residue's atoms
        for (size_t j = other->m_count; j > 0; j--) {
            // get the other atom
            PROT::Atom * other_atom = &other->m_atoms[j-1];
            // if the other atom is a hydrogen, skip it
            // if (other_atom->m_element == "H" or other_atom->is_backbone_atom()) {
            if (other_atom->m_element == "H") {
                continue;
            }
            // check for a clash
            if (atom->clash(other_atom)) {
                return true;
            }
        }
    }
    return false;
}
