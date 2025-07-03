/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * score method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Calculate a score for the Residue that tells how "complete" it is
double PROT::Residue::score () const {
    // If the Residue is missing or has no Atoms, the score is 0
    if ((!m_present) || (m_count == 0)) {return 0.0;}
    // If the Residue is present and has no missing Atoms, the score is 1
    else if (m_missing_atoms == 0) {return 1.0;}
    // Otherwise, calculate the score as follows
    return ((double) m_count / ((double) (m_count + m_missing_atoms)));
}
