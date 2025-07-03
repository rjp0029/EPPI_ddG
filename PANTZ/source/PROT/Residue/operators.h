/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * [] operators of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Access to the Atoms of the Residue
PROT::Atom * PROT::Residue::get_atom (const size_t i) {
    // If there are no atoms, raise an error
    if (m_count == 0) {
        string error = "It is not possible to access an Atom in a Residue "
                       "when the Residue is empty.\n";
        throw PANTZ_error (error);}
    // Also raise an error if the index is too large
    else if (i >= m_count) {
        stringstream c1; c1 << m_count;
        stringstream c2; c2 << i;
        string error = "This Residue contains " + c1.str() + " Atoms. "
                     + c2.str() + " is not a valid index.\n";
        for(size_t j=0; j<m_count; ++j) {error += m_atoms[j].str();}
        throw PANTZ_error (error);}
    // Return a pointer to the requested atom
    return &(m_atoms[i]);
}

// Access Atoms using a string
PROT::Atom * PROT::Residue::get_atom (const string& label) {
    // If there are no atoms, raise an error
    if (m_count == 0) {
        string error = "It is not possible to access an Atom in a Residue "
                       "when the Residue is empty.\n";
        throw PANTZ_error (error);}
    // Find the atom
    for(size_t i=0; i<m_count; ++i) {
        if (m_atoms[i].m_name == label) {return &(m_atoms[i]);}}
    // If there was no such atom
    string error = "This Residue does not contain a " + label + " Atom.\n";
    for (size_t i=0; i<m_count; ++i) {error += m_atoms[i].str();}
    throw PANTZ_error (error);
    return 0;
}

// Access to AtomPtrs instead of pointers to Atoms
PROT::AtomPtr PROT::Residue::operator[] (const size_t i) {
    AtomPtr ptr = get_atom(i);
    return ptr;
}

PROT::AtomPtr PROT::Residue::operator[] (const string& label) {
    AtomPtr ptr = get_atom(label);
    return ptr;
}

bool PROT::Residue::operator<(const PROT::Residue& other) const {
    // compare based on number of atoms
    if (m_count < other.m_count) {
        return true;
    }
    // otherwise return false
    return false;
}
