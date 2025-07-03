/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the allocate methods of the class. */

// Error check the inclusion chain for the file. Note that this file is included
// by Atom.h, not Matrix.h, because it requires use of Atom methods. Atom.h
// checks that it is loaded after Matrix.h, so that is fine.
#ifndef Atom_Loading_Status
#error The allocate_atom.h Matrix file must be included by Atom.h
#endif

// Implement the methods

// Allocate a Matrix from an Atom's coordinates
void PROT::Matrix::allocate (const Atom * atom) {
    // Delete existing information
    clean_up();
    // there is 1 row
    m_rows = 1;
    // And 3 columns
    m_columns = AtomCoordinates;
    // The number of values in the matrix
    size_t n = m_rows * m_columns;
    // Allocate the memory
    m_values = new coor [n];
    // Assign the atom's coordinates as the values
    for (size_t i=0; i<n; ++i) {m_values[i] = atom->operator[](i);}
}

// Allocate a Matrix from a vector of Atoms
void PROT::Matrix::allocate(const vector<Atom *>& atoms) {
    // Delete existing information
    clean_up ();
    // Throw an error if the vector is empty
    if (atoms.size() == 0) {
        string error = "A Matrix cannot be allocated from an empty vector of "
                       "Atoms.\n";
        throw PANTZ_error (error);}
    // Set the number of rows
    m_rows = atoms.size();
    // The number of columns
    m_columns = AtomCoordinates;
    // The total number of values
    size_t n = m_rows * m_columns;
    // Allocate the memory
    m_values = new coor [n];
    // A counter of which value is being stored
    size_t p = 0;
    // Go through the atoms and their coordinates
    for(size_t i=0; i<atoms.size(); ++i) {
        for (size_t j=0; j<AtomCoordinates; ++j) {
            m_values[p] = atoms[i]->operator[](j);
            ++p;}}
}

// Allocate a Matrix from a single AtomPtr
void PROT::Matrix::allocate (AtomPtr& atom) {allocate(atom.pointer());}

// Allocate a Matrix from a vector of AtomPtrs
void PROT::Matrix::allocate (vector<AtomPtr>& atoms) {
    // Create a vector of atom pointers instead
    vector<Atom *> ptrs; ptrs.reserve(atoms.size());
    for(size_t i=0; i<atoms.size(); ++i) {
        ptrs.push_back(atoms[i].pointer());}
    // Now allocate using those
    allocate(ptrs);
}
