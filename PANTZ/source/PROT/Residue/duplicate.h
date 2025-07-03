/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * duplicate method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Create another Residue that has the same contents as this Residue
PROT::Residue PROT::Residue::duplicate () const {
    // Create the other Residue
    Residue other;
    // Copy the contents of this Residue to that Residue
    other.m_count = m_count;
    other.m_name = m_name;
    other.m_number = m_number;
    other.m_internal = m_internal;
    other.m_insertion = m_insertion;
    other.m_protein = m_protein;
    other.m_present = m_present;
    other.m_missing_atoms = m_missing_atoms;
    other.m_phi = m_phi;
    other.m_psi = m_psi;
    other.m_omega = m_omega;
    // If there are Atoms, allocate them
    if (m_count > 0) {
        other.m_atoms = new PROT::Atom [m_count];
        for(size_t i=0; i<m_count; ++i) {
            other.m_atoms[i] = m_atoms[i];}}
    // Return the other residue
    return other;
}
