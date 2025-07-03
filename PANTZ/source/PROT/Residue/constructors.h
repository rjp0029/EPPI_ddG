/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * construction methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// construct a residue from a vector of atoms
PROT::Residue::Residue (const vector<PROT::Atom *>& atoms) {
    // Initialize the residue
    initialize();
    // do a complete load of the Residue's atoms
    load(atoms, false, true);
}

// Construct a residue from a vector of AtomPtrs
PROT::Residue::Residue (vector<AtomPtr>& atoms) {
    // Initialize the residue
    initialize ();
    // Do a complete load of the Residue's atoms
    load (atoms, false, true);
}

// Construct a Residue from a string name and the protein's character. This
// function should really only be used in the context of a PDB file
PROT::Residue::Residue (const string& input, const char L) {
    // Initialize the Residue
    initialize();
    // Set the Residue's name
    m_name = input;
    // Set the Protein's name
    m_protein = L;
    // The remaining attributes can be left as their initialized values
}
