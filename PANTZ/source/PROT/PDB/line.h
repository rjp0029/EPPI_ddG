/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the line method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Exterior access to a line in the PDB file
string PROT::PDB::line (const size_t i) const {
    if (i >= m_lines.size()) {
        stringstream c1; c1 << m_lines[i].size();
        stringstream c2; c2 << i;
        string error = c1.str() + " lines were loaded from PDB file "
                     + m_name + ". An index of " + c2.str() + " is not "
                       "acceptable to access one of them.\n";
        throw PANTZ_error (error);}
    return m_lines[i];
}
