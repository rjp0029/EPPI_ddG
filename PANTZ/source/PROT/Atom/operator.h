/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the [] operator methods of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Provide access to the Atom's coordinates using an index. Because size_t is an
// unsigned integer, the indexing must be between 0 and AtomCoordinates
PROT::coor PROT::Atom::operator[] (const size_t i) const {
    // If the index is not valid, throw an error
    if (i >= AtomCoordinates) {
        stringstream c1; c1 << AtomCoordinates;
        stringstream c2; c2 << i;
        string error = "Atoms have " + c1.str() + " coordinates.\n"
                     + c2.str() + " is not an acceptable value.\n";
        throw PANTZ_error (error);}
    // Return the coordinate
    return m_coors[i];
}

// Provide access to the Atom's coordinates using a character
PROT::coor PROT::Atom::operator[] (const char L) const {
    // Permit access to the x, y, and z coordinates
    if ((L == 'x') || (L == 'X')) {return m_coors[0];}
    else if ((L == 'y') || (L == 'Y')) {return m_coors[1];}
    else if ((L == 'z') || (L == 'Z')) {return m_coors[2];}
    // throw an error
    string error = "'";
    error += L;
    error += "' is not a recognized Atom coordinate.\n";
    throw PANTZ_error(error);
    // Return a generic value so the function compiles correctly
    return 0.0;
}


