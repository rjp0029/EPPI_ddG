/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the find_neighbors functions of the NeighborBox class.
 * This is the class that actually finds the nearby atoms to the specified point
 * */

// Make sure it is being included as expected
#ifndef NeighborBox_Loading_Status
#error The NeighborBox find file must be included by NeighborBox.h
#endif

// Implement the method that fills a vector of Atoms
void PROT::NeighborBox::find_neighbors (vector<Atom *>& options,
     const coor X, const coor Y, const coor Z, const coor D) {
    // Clear the options vector
    options.clear();
    // Calculate the minimum and maximum permissible X, Y, and Z coordinates
    coor xmin = X - D;
    coor xmax = X + D;
    coor ymin = Y - D;
    coor ymax = Y + D;
    coor zmin = Z - D;
    coor zmax = Z + D;
    // Calculate the bins in which to search
    size_t XL = calculate_bin(xmin, 'X');
    size_t XU = calculate_bin(xmax, 'X');
    size_t YL = calculate_bin(ymin, 'Y');
    size_t YU = calculate_bin(ymax, 'Y');
    size_t ZL = calculate_bin(zmin, 'Z');
    size_t ZU = calculate_bin(zmax, 'Z');
    // Search for atoms in the identified bins
    for(size_t xb=XL; xb<=XU; ++xb) {
        for(size_t yb=YL; yb<=YU; ++yb) {
            for(size_t zb=ZL; zb<=ZU; ++zb) {
                // Loop through the atoms in the bin
                for(size_t i=0; i<m_atom_pointers[xb][yb][zb].size(); ++i) {
                    // Get the atom pointer
                    Atom * atom = m_atom_pointers[xb][yb][zb][i];
                    // Get it's x coordinate
                    size_t c = 0; coor x = (*atom)[c];
                    // If it is outside of either bound, skip it
                    if ((x < xmin) || (x > xmax)) {continue;}
                    // Do the same for y and z
                    c = 1; coor y = (*atom)[c];
                    if ((y < ymin) || (y > ymax)) {continue;}
                    c = 2; coor z = (*atom)[c];
                    if ((z < zmin) || (z > zmax)) {continue;}
                    // Since the atom's coordinates are in the required range,
                    // store it
                    options.push_back(atom);}}}}
    // End the function
}

// Do the same thing with filling a vector of Atom Pointers instead
void PROT::NeighborBox::find_neighbors(vector<AtomPtr>& options,
     const coor X, const coor Y, const coor Z, const coor D) {
    // Make a vector of pointers to atoms
    vector<Atom *> atoms;
    // Find neighbors using that
    find_neighbors(atoms, X, Y, Z, D);
    // Use that to fill in the options
    options.clear(); options.reserve(atoms.size());
    for(size_t i=0; i<atoms.size(); ++i) {
        options.push_back(AtomPtr(atoms[i]));}
}
