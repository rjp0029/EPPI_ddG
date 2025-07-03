/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the initialize method of the Atom class. */

// Confirm that the Atom class has been and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Assign default values to the Atom's attributes
void PROT::Atom::initialize () {
    // The atom's coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = 0.0;}
    // String attributes
    m_type = "ATOM";
    m_name = "NONE";
    m_residue = "NON";
    m_element = "";
    m_charge = "";
    // Integer attributes
    m_number = 1;
    m_residue_number = 1;
    // Decimal attributes
    m_occupancy = 1.0;
    m_temperature = 0.0;
    // Char attributes
    m_alt = ' ';
    m_insertion = ' ';
    m_protein = ' ';
}

