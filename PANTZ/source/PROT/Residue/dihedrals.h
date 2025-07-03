/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * phi, psi, and omega methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Access to the Residue's dihedral angles
PROT::coor PROT::Residue::phi () const {
    // Throw an error if the angle was never assigned
    if (m_phi < -999) {
        string error = "No phi angle available.\n";
        throw PANTZ_error (error);}
    return m_phi;
}

PROT::coor PROT::Residue::psi () const {
    if (m_psi < -999) {
        string error = "No psi angle available.\n";
        throw PANTZ_error (error);}
    return m_psi;
}

PROT::coor PROT::Residue::omega () const {
    if (m_omega < -999) {
        string error = "No omega angle available.\n";
        throw PANTZ_error (error);}
    return m_omega;
}
