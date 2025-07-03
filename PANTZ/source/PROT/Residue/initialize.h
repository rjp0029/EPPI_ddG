/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * initialize method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The initialization function of the Residue class
void PROT::Residue::initialize () {
    // Assign default values to all variables
    m_atoms = 0;
    m_count = 0;
    m_name = "N/A";
    m_number = 0;
    m_internal = 0;
    m_insertion = ' ';
    m_protein = ' ';
    m_present = true;
    m_missing_atoms = 0;
    m_phi = -1000.0;
    m_psi = -1000.0;
    m_omega = -1000.0;
}
