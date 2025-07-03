/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the method to determine if two atoms are sterically clashing based on
 * 90% vdw overlap. 
 *
*/

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// determine if this atoms is clashing with another
bool PROT::Atom::clash(Atom* other, float probe_radius) {
    // get the vdw radii
    float vdw1 = this->lj_sigma() + probe_radius;
    float vdw2 = other->lj_sigma() + probe_radius;
    // get the distance between the atoms
    float distance = this->distance(*other);
    // check for a clash
    if (distance < 0.9*(vdw1 + vdw2)) {
        return true;
    }
    return false;
}

