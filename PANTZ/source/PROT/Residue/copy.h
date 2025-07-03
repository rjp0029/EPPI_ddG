/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * pointer-based copy method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Copy the information from another Residue into this one
void PROT::Residue::copy (const Residue * other) {
    // Clean up any existing information
    clean_up();
    // Get the number of Atoms in the other residue
    m_count = other->m_count;
    // If there are Atoms, allocate and copy them
    if (m_count > 0) {
        m_atoms = new PROT::Atom [m_count];
        for(size_t i=0; i<m_count; ++i) {m_atoms[i] = other->m_atoms[i];}}
    // Copy the other attributes
    m_name = other->m_name;
    m_number = other->m_number;
    m_internal = other->m_internal;
    m_insertion = other->m_insertion;
    m_protein = other->m_protein;
    m_present = other->m_present;
    m_missing_atoms = other->m_missing_atoms;
    m_phi = other->m_phi;
    m_psi = other->m_psi;
    m_omega = other->m_omega;
}
