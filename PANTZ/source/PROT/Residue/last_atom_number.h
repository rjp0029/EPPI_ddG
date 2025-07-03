/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to get the number of the last atom in the Residue. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Get the number of the last Atom in the Residue
long PROT::Residue::last_atom_number () const {
    // If there are atoms
    if (m_count > 0) {return m_atoms[m_count-1].m_number;}
    // Otherwise return 0
    return 0;
}
