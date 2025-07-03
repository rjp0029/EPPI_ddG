/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the construction methods of the NeighborBox class. */

// Make sure it is being included as expected
#ifndef NeighborBox_Loading_Status
#error The NeighborBox construction file must be included by NeighborBox.h
#endif

// Implement the private construction method that all of the other methods use
void PROT::NeighborBox::construct (vector<Atom *>& atoms, const coor BL,
                                   const bool considerHydrogen) {
    // Initialize the class
    initialize ();
    // Store the bin side length
    m_bin_side = BL;
    // If the value isn't reasonable, throw an error
    if (m_bin_side <= 0) {
        string error = "The NeighborBox Bin Length must be a positive number.\n";
        throw PANTZ_error(error);}
    // Loop through the Atoms
    for(size_t i=0; i<atoms.size(); ++i) {
        // Determine if this atom is worth considering
        if ((considerHydrogen) || (!atoms[i]->is_hydrogen())) {
            // Get the X, Y, and Z coordinates of the atom
            size_t c = 0; coor x = atoms[i]->operator[](c);
            c = 1; coor y = atoms[i]->operator[](c);
            c = 2; coor z = atoms[i]->operator[](c);
            // If this is the first atom being considered, store the coordinates
            // in the bounds
            if (m_size == 0) {
                m_XL = x; m_XU = x;
                m_YL = y; m_YU = y;
                m_ZL = z; m_ZU = z;}
            // Otherwise update them if they exceed the bounds
            else {
                if (x < m_XL) {m_XL = x;}
                if (x > m_XU) {m_XU = x;}
                if (y < m_YL) {m_YL = y;}
                if (y > m_YU) {m_YU = y;}
                if (z < m_ZL) {m_ZL = z;}
                if (z > m_ZU) {m_ZU = z;}}
            // Increment the atom counter
            ++m_size;}}
    // If there are no atoms to consider, throw an error
    if (m_size == 0) {
        string error = "NeighborBox constructed to store no atoms.\n";
        throw PANTZ_error (error);}
    // Calculate the number of x, y, and z bins
    m_Xbins = 1 + (size_t) ((m_XU - m_XL) / m_bin_side);
    m_Ybins = 1 + (size_t) ((m_YU - m_YL) / m_bin_side);
    m_Zbins = 1 + (size_t) ((m_ZU - m_ZL) / m_bin_side);
    // Allocate the memory for the bins
    m_atom_pointers.resize(m_Xbins);
    for(size_t i=0; i<m_Xbins; ++i) {
        m_atom_pointers[i].resize(m_Ybins);
        for(size_t j=0; j<m_Ybins; ++j) {
            m_atom_pointers[i][j].resize(m_Zbins);}}
    // Loop through the atoms and put them in bins
    for(size_t i=0; i<atoms.size(); ++i) {
        if ((considerHydrogen) || (!atoms[i]->is_hydrogen())) {
            size_t c = 0; size_t x = calculate_bin(atoms[i]->operator[](c), 'X');
            c = 1; size_t y = calculate_bin(atoms[i]->operator[](c), 'Y');
            c = 2; size_t z = calculate_bin(atoms[i]->operator[](c), 'Z');
            m_atom_pointers[x][y][z].push_back(atoms[i]);}}
    // End the function
}

// Now implement the six public constructors of the class
PROT::NeighborBox::NeighborBox (vector<Atom>& atoms, const coor BL, 
      const bool considerHydrogen = true) {
    // Make a vector of pointers to the atoms
    vector<Atom *> use; use.reserve(atoms.size());
    for(size_t i=0; i<atoms.size(); ++i) {use.push_back(&atoms[i]);}
    // Call the construct function
    construct(use, BL, considerHydrogen);
}

PROT::NeighborBox::NeighborBox (vector<AtomPtr>& atoms, const coor BL,
      const bool considerHydrogen = true) {
    vector<Atom *> use; use.reserve(atoms.size());
    for(size_t i=0; i<atoms.size(); ++i) {use.push_back(atoms[i].pointer());}
    construct(use, BL, considerHydrogen);
}

PROT::NeighborBox::NeighborBox (vector<Residue>& residues, const coor BL,
      const bool considerHydrogen = true) {
    // Make the vector of atom pointers
    vector<Atom *> use;
    // Loop through the residues
    for(size_t i=0; i<residues.size(); ++i) {
        // Loop through the Residue's atoms
        for(size_t j=0; j<residues[i].size(); ++j) {
            // Store the atom
            use.push_back(residues[i].get_atom(j));}}
    // Call the constructor
    construct(use, BL, considerHydrogen);
}

PROT::NeighborBox::NeighborBox (vector<ResiduePtr>& residues, const coor BL,
      const bool considerHydrogen = true) {
    vector<Atom *> use;
    for(size_t i=0; i<residues.size(); ++i) {
        for(size_t j=0; j<residues[i].size(); ++j) {
            use.push_back(residues[i].get_atom(j));}}
    construct(use, BL, considerHydrogen);
}

PROT::NeighborBox::NeighborBox (vector<Protein>& proteins, const coor BL,
      const bool considerHydrogen = true) {
    // Make the vector of atoms
    vector<Atom *> use;
    // Loop through the proteins
    for(size_t i=0; i<proteins.size(); ++i) {
        // Loop through it's residues
        for(size_t j=0; j<proteins[i].size(); ++j) {
            // Loop through the atoms
            for(size_t k=0; k<proteins[i][j].size(); ++k) {
                use.push_back(proteins[i][j].get_atom(k));}}}
    // Call the constructor
    construct(use, BL, considerHydrogen);
}

PROT::NeighborBox::NeighborBox (vector<Protein>& prots1, 
      vector<Protein>& prots2, const coor BL, 
      const bool considerHydrogen = true) {
    vector<Atom *> use;
    for(size_t i=0; i<prots1.size(); ++i) {
        for(size_t j=0; j<prots1[i].size(); ++j) {
            for(size_t k=0; k<prots1[i][j].size(); ++k) {
                use.push_back(prots1[i][j].get_atom(k));}}}
    for(size_t i=0; i<prots2.size(); ++i) {
        for(size_t j=0; j<prots2[i].size(); ++j) {
            for(size_t k=0; k<prots2[i][j].size(); ++k) {
                use.push_back(prots2[i][j].get_atom(k));}}}
    construct(use, BL, considerHydrogen);
}
