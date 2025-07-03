/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * rmsd methods of the residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The function to calculate the rmsd between this residue and another
float PROT::Residue::rmsd(PROT::Residue* other) {
    // Get the atoms of this residue
    vector<PROT::Atom*> atoms;
    // Get the atoms of the other residue
    vector<PROT::Atom*> other_atoms;
    // go through each set of atoms and add them to the vectors enruing that they are the same order (by name)
    for (size_t i = 0; i < this->size(); i++) {
        for (size_t j = 0; j < other->size(); j++) {
            if (this->m_atoms[i].name() == other->m_atoms[j].name()) {
                atoms.push_back(&this->m_atoms[i]);
                other_atoms.push_back(&other->m_atoms[j]);
            }
        }
    }
    // Make sure the two residues have the same number of atoms
    if (atoms.size() != other_atoms.size()) {
        string error = "The two residues do not have common number of atoms.";
        throw PANTZ_error(error);
    }
    // Calculate the rmsd
    float rmsd = 0;
    for (size_t i = 0; i < atoms.size(); i++) {
        rmsd += pow(atoms[i]->distance(*other_atoms[i]), 2);
    }
    rmsd = sqrt(rmsd / atoms.size());
    return rmsd;
}
