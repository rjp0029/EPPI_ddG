/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * center method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Move a Residue so that its center of mass is at the origin
PROT::Matrix PROT::Residue::center (const string how = "heavy") {
    // Allocate a 1 x AtomCoordinates matrix
    Matrix output (1, AtomCoordinates);
    // If there are atoms
    if (m_count > 0) {
        // Select the specified atoms
        vector<Atom *> chosen; select_atoms(chosen, how);
        // If none were chosen, throw an error
        if (chosen.size() == 0) {
            string error = "This Residue could not be centered using the "
                         + how + " Atom selection method:\n" + str();
            throw PANTZ_error (error);}
        // Loop through the coordinates of Atoms
        for(size_t i=0; i<AtomCoordinates; ++i) {
            // Add up the values of the Atoms for this dimension
            coor n = 0.0;
            for(size_t j=0; j<chosen.size(); ++j) {
                n += chosen[j]->m_coors[i];}
            // Calculate the average value
            n /= m_count;
            // Store it in the matrix
            output.set(0, i, n);}
        // Use the private move function of the class to subtract those
        // coordinates from every atom in the Residue
        private_move(&output, '-');}
    return output;
}

// Same thing from a ResiduePtr
PROT::Matrix PROT::ResiduePtr::center (const string how = "heavy") {
    check ();
    return m_ptr->center(how);
}
