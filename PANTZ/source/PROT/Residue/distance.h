/* Created by the clay at Auburn University
 *
 * This file is intended to be included by the Residue.h header file. It contains
 * the method to get the distance to antother atom. */

// Confirm that the Residue class has been declared and is actively being loaded
#ifndef Residue_Loading_Status
#error Residue methods must be included from the Residue.h header file
#endif 

// get the distance between CA atoms of two residues
float PROT::Residue::distance (Residue& other) {
    // get the distance from this CA atom to the other CA atom
    float distance = get_atom("CA")->distance(*other.get_atom("CA"));
    return distance;
}

// get the minimum distance between sets of atoms on two residues
float PROT::Residue::min_distance (Residue& other) {
    // set initial to float max
    float min_distance = std::numeric_limits<float>::max();
    // loop through the atoms in this residue
    for (size_t i = 0; i < size(); i++) {
        // get the atom
        Atom* atom1 = get_atom(i);
        // loop through the atoms in the other residue
        for (size_t j = 0; j < other.size(); j++) {
            // get the atom
            Atom* atom2 = other.get_atom(j);
            // get the distance between the atoms
            float distance = atom1->distance(*atom2);
            // if the distance is less than the minimum distance, set the minimum distance
            if (distance < min_distance) {
                min_distance = distance;
            }
        }
    }
    return min_distance;
}