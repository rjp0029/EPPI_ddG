/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * position method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Position a Residue in a standard way to facilitate rotamer positioning
void PROT::Residue::position (PROT::Matrix& how, const bool forward = true) {
    // Do one set of operations if the Residue is being assigned the initial
    // position
    if (forward) {
        // This function only works on Amino Acid residues (permit Histidine
        // variations, too)
        if (!is_amino_acid()) {
            // Throw an error
            string error = "The Residue position function only works on amino "
                           "acids, not " + m_name + ".\n";
            throw PANTZ_error (error);}
        // The three Atoms that are used in the positioning. Store them here
        vector<Atom *> atoms; atoms.reserve(3); select_atoms(atoms, "rotamer");
        // Allocate the how matrix to be the proper dimensions
        how.allocate(1, AtomCoordinates + 3);
        // Create a matrix using the coordinates of the first atom
        Matrix first (atoms[0]);
        // Move the Residue so that first atom is at the origin
        private_move(&first, '-');
        // Store the coordinates in the how matrix
        for(size_t i=0; i<AtomCoordinates; ++i) {how.set(0, i, first(0, i));}
        // Calculate the angle of rotation to rotate the second Atom into the
        // XZ plane
        coor angle1 = fabs(acos(atoms[1]->m_coors[0] /
                           sqrt(pow(atoms[1]->m_coors[0], 2) 
                              + pow(atoms[1]->m_coors[1], 2))));
        if (atoms[1]->m_coors[1] > 0) {angle1 = -angle1;}
        // A unit vector for the Z axis
        coor Z [3] = {0.0, 0.0, 1.0};
        // Create a rotation matrix
        Matrix second (angle1, Z);
        // Rotate the Residue
        private_rotate(&second);
        // Store the angle in the how matrix
        how.set(0, AtomCoordinates, angle1);
        // Rotate the second Atom onto the Z-axis by rotating around the
        // y-axis
        coor angle2 = fabs(acos(atoms[1]->m_coors[2] /
                           sqrt(pow(atoms[1]->m_coors[2], 2) 
                              + pow(atoms[1]->m_coors[0], 2))));
        if (atoms[1]->m_coors[0] > 0) {angle2 = -angle2;}
        coor Y [3] = {0.0, 1.0, 0.0};
        Matrix third (angle2, Y);
        private_rotate(&third);
        how.set(0, AtomCoordinates+1, angle2);
        // Rotate the third atom into the XZ plane
        coor angle3 = fabs(acos(atoms[2]->m_coors[0] / 
                           sqrt(pow(atoms[2]->m_coors[0], 2)
                              + pow(atoms[2]->m_coors[1], 2))));
        if (atoms[2]->m_coors[1] > 0) {angle3 = -angle3;}
        Matrix fourth (angle3, Z);
        private_rotate(&fourth);
        how.set(0, AtomCoordinates+2, angle3);
        // The residue is now positioned
    }
    // If this is not the initial positioning, reverse the calculations
    else {
        coor Z [3] = {0.0, 0.0, 1.0}; coor Y [3] = {0.0, 1.0, 0.0};
        Matrix fourth (-how(0, AtomCoordinates+2), Z);
        private_rotate(&fourth);
        Matrix third (-how(0, AtomCoordinates+1), Y);
        private_rotate(&third);
        Matrix second (-how(0, AtomCoordinates), Z);
        private_rotate(&second);
        Matrix first (1, AtomCoordinates);
        for(size_t i=0; i<AtomCoordinates; ++i) {first.set(0, i, how(0, i));}
        private_move(&first, '+');}
}

// The same function in the ResiduePtr class
void PROT::ResiduePtr::position (Matrix& matrix, const bool forward = true) {
    check ();
    m_ptr->position(matrix, forward);
}
