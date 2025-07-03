/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * str method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Generate a string of the Residue's information properly formatted for
// internal or external use
string PROT::Residue::str (bool internal = false) {
    // Store the output string here
    string output;
    // If there are atoms
    if (m_count > 0) {
        // Allocate an appropriate amount of space in the string
        output.reserve(AtomStringLength * m_count);
        // Loop through the Atoms
        for(size_t i=0; i<m_count; ++i) {
            // If internal numbering should be used
            if (internal) {
                m_atoms[i].m_residue_number = m_internal;
                m_atoms[i].m_insertion = ' ';}
            else {
                m_atoms[i].m_residue_number = m_number;
                m_atoms[i].m_insertion = m_insertion;}
            // Add the Atom's string to the output
            output.append(m_atoms[i].str());}}
    return output;
}

// The ResiduePtr version
string PROT::ResiduePtr::str (bool internal = false) {
    check();
    return m_ptr->str(internal);
}
