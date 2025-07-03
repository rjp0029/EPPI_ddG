/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * clean_up method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Delete any dynamically allocated memory
void PROT::Residue::clean_up () {
    if (m_atoms != 0) {delete[] m_atoms; m_atoms = 0;}
}
