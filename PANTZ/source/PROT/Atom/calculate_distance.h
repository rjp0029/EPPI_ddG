/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the calculate_distance methods of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Calculate the distance between 2 atoms
PROT::coor PROT::Atom::calculate_distance (const Atom * other, 
                                           const bool squared = false) const {
    // The distance will be stored here
    coor value = 0.0;
    // Loop through the coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Add the square of the difference between the coordinates
        value += pow(m_coors[i] - other->m_coors[i], 2);}
    // If appropriate, return just the squared sum
    if (squared) {return value;}
    // Calculate the square root of the sum
    return sqrt(value);
}

// Do the same calculations, just using an Atom passed by reference instead
PROT::coor PROT::Atom::calculate_distance (const Atom& other, 
                                           const bool squared = false) const {
    return calculate_distance(&other, squared);
}

// Calculate distances using AtomPtrs 
PROT::coor PROT::AtomPtr::calculate_distance (const Atom * other,
                                              const bool squared = false) const {
    check ();
    return m_ptr->calculate_distance (other, squared);
}

PROT::coor PROT::AtomPtr::calculate_distance (const Atom& other,
                                              const bool squared = false) const {
    check ();
    return m_ptr->calculate_distance (other, squared);
}

PROT::coor PROT::AtomPtr::calculate_distance (const AtomPtr& other,
                                              const bool squared = false) const {
    check ();
    other.check();
    return m_ptr->calculate_distance(other.m_ptr, squared);
}
