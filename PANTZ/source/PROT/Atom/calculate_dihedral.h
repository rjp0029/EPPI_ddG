/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the calculate_dihedral method of the PROT namespace. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Calculate the dihedral angle between 4 atoms
PROT::coor PROT::calculate_dihedral (const Atom * atom1, const Atom * atom2,
                                     const Atom * atom3, const Atom * atom4) {
    // Make vectors of the differences between certain atoms
    Matrix f = Matrix(atom1) - Matrix(atom2); f.make_unit_vector();
    Matrix g = Matrix(atom2) - Matrix(atom3); g.make_unit_vector();
    // This is supposed to be 4 - 3
    Matrix h = Matrix(atom4) - Matrix(atom3); h.make_unit_vector();
    // Calculate cross products between the vectors
    Matrix a = f.crossProduct(g); a.make_unit_vector();
    Matrix b = h.crossProduct(g); b.make_unit_vector(); b = b.transpose();
    // Calculate the dot product of a and b
    Matrix c = a.dotProduct(&b);
    // Calculate the angle
    coor angle = (180.0/M_PI) * acos(c(0,0));
    // Calculate appropriate values for checking if the sign on the angle
    // needs to be changed
    Matrix check1 = a.crossProduct(b); check1 = check1.transpose();
    Matrix check2 = g.dotProduct(&check1);
    if (check2(0,0) > 0) {angle = -angle;}
    return angle;
}

// Do the same, using 4 atoms passed by reference
PROT::coor PROT::calculate_dihedral (const Atom& atom1, const Atom& atom2,
                                     const Atom& atom3, const Atom& atom4) {
    return calculate_dihedral(&atom1, &atom2, &atom3, &atom4);
}

// Do the same using 4 AtomPtr objects
PROT::coor PROT::calculate_dihedral (AtomPtr& atom1, AtomPtr& atom2,
                                     AtomPtr& atom3, AtomPtr& atom4){
    // Get the 4 pointers, which also checks to make sure they aren't null
    const Atom * ptr1 = atom1.pointer();
    const Atom * ptr2 = atom2.pointer();
    const Atom * ptr3 = atom3.pointer();
    const Atom * ptr4 = atom4.pointer();
    // Calculate the angle
    return calculate_dihedral (ptr1, ptr2, ptr3, ptr4);
}    
