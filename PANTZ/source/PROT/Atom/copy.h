/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the pointer-based copy method of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Copy the information from another Atom into this atom
void PROT::Atom::copy (const Atom * other) {
    // Copy all of the attributes
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = other->m_coors[i];}
    m_type = other->m_type;
    m_name = other->m_name;
    m_residue = other->m_residue;
    m_element = other->m_element;
    m_charge = other->m_charge;
    m_number = other->m_number;
    m_residue_number = other->m_residue_number;
    m_occupancy = other->m_occupancy;
    m_temperature = other->m_temperature;
    m_alt = other->m_alt;
    m_insertion = other->m_insertion;
    m_protein = other->m_protein;
}

