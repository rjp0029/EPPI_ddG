/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the structure method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Access to a Structure in the PDB file
PROT::Structure * PROT::PDB::structure (const size_t i) {
    if (i >= m_structures.size()) {
        stringstream c1; c1 << m_structures[i].proteins();
        stringstream c2; c2 << i;
        string error = c1.str() + " Structures were identified in PDB file "
                     + m_name + ". An index of " + c2.str() + " is not "
                       "accepatble to access one of them.\n";
        throw PANTZ_error (error);}
    return &(m_structures[i]);
}
