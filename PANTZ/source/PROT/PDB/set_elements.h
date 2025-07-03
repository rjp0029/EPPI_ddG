/* Created by the clay at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the method to set the elements if they are not already set. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// set the m_elements of the atoms in the PDB file
void PROT::PDB::set_elements () {
    // Set the elements given the name, after updating atoms after rosetta
    for (size_t i = 0; i < m_proteins.size(); i++) {
        m_proteins[i].update_atoms_after_Rosetta();
        for (size_t j = 0; j < m_proteins[i].size(); j++) {
            for (size_t k = 0; k < m_proteins[i](j, ' ', true)->size(); k++) {
                char element = m_proteins[i](j, ' ', true)->get_atom(k)->determine_element();
                m_proteins[i](j, ' ', true)->get_atom(k)->m_element = element;
            }
        }
    }
}
