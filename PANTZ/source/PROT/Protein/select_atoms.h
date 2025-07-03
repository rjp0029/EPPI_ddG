/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the atom selection methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Select atoms and store them as pointers
void PROT::Protein::select_atoms (vector<Atom *>& atoms, 
                                  const string how = "all") {
    // Make sure the protein is not empty
    if (m_count == 0) {
        string error = "It is not possible to select atoms from an empty "
                       "Protein.\n";
        throw PANTZ_error (error);}
    // Go through the residues
    for(size_t i=0; i<m_count; ++i) {
        // Store atoms from this residue
        m_residues[i].select_atoms(atoms, how);}
}

// Store them as AtomPtrs
void PROT::Protein::select_atoms (vector<AtomPtr>& atoms,
                                  const string how = "all") {
    // Make a vector of atom pointers
    vector<Atom *> chosen;
    // Select the atoms 
    select_atoms(chosen, how);
    // Store them in atoms
    for(size_t i=0; i<chosen.size(); ++i) {
        atoms.push_back(chosen[i]);}
}