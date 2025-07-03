/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the implementation of
 * the residue allocation methods of the class. */

// Error check the inclusion chain for the file. Note that this file is included
// by Residue.h, not Matrix.h, because it requires use of Residue methods. 
// Residue.h checks that it is loaded after Matrix.h, so that is fine.
#ifndef Residue_Loading_Status
#error The allocate_residue.h Matrix file must be included by Atom.h
#endif

// Implement the methods

// Allocate a Matrix from a Residue
void PROT::Matrix::allocate (Residue * res, const string how = "heavy") {
    // Get the atoms to use from the residue
    vector<Atom *> atoms; res->select_atoms(atoms, how);
    // Allocate the matrix from those atoms
    try {allocate(atoms);}
    catch (PANTZ_error& e) {
        string error = "This error occurred because the " + how + " atom "
                       "selection method found no atoms for this Residue:\n"
                     + res->str();
        throw PANTZ_error (e, error);}
}

// Do the same with a referenced matrix
void PROT::Matrix::allocate (Residue& res, const string how = "heavy") {
    allocate(&res, how);
}

// And a Residue Pointer
void PROT::Matrix::allocate (ResiduePtr& res, const string how = "heavy") {
    allocate(res.pointer(), how);
}

// Matrix constructors from Residues
PROT::Matrix::Matrix (Residue * res, const string how = "heavy") {
    initialize(); allocate(res, how);
}

PROT::Matrix::Matrix (Residue& res, const string how = "heavy") {
    initialize(); allocate(res, how);
}

PROT::Matrix::Matrix (ResiduePtr& res, const string how = "heavy") {
    initialize(); allocate(res, how);
}
