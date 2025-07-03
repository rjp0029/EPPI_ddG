/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * atom renumbering method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Renumber the Atoms in the Residue
long PROT::Residue::renumber_atoms (long n) {
    // If there are Atoms
    if (m_count > 0) {
        // Loop through them
        for(size_t i=0; i<m_count; ++i) {
            // modify the atom's number
            m_atoms[i].m_number = n;
            // Increment n
            ++n;}}
    return n;
}
