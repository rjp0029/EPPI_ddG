/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * check_atoms method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Determine whether or not a vector of Atom pointers are acceptable for use
// in a Residue
void PROT::Residue::check_atoms (const vector<PROT::Atom *>& atoms) const {
    // If there are no Atoms, throw an error
    if (atoms.size() == 0) {
        string error = "A Residue cannot be assembled from an empty vector "
                       "of Atoms.\n";
        throw PANTZ_error (error);}
    // If there is only a single Atom, there cannot be naming / numbering
    // issues with it
    else if (atoms.size() == 1) {return;}
    // Make sure that every Atom has the same residue and molecule information
    for (size_t i=1; i<atoms.size(); ++i) {
        if (((atoms[0]->m_residue != atoms[i]->m_residue) || 
             (atoms[0]->m_residue_number != atoms[i]->m_residue_number)) ||
            ((atoms[0]->m_insertion != atoms[i]->m_insertion) ||
             (atoms[0]->m_protein != atoms[i]->m_protein))) {
            string error = "All Atoms in a Residue must have the same "
                           "residue and protein labelling information. These "
                           "do not:\n";
            for(size_t j=0; j<atoms.size(); ++j) {
                error += atoms[j]->str();}
            throw PANTZ_error (error);}}
    // Confirm that each ATOM entry has a unique name
    for(size_t i=0; i<atoms.size()-1; ++i) {
        // Skip HETATM entries
        if (atoms[i]->m_type == "HETATM") {continue;}
        // Loop through all subsequent atoms
        for(size_t j=i+1; j<atoms.size(); ++j) {
            // Skip HETATM entries
            if (atoms[j]->m_type == "HETATM") {continue;}
            // Throw an error if the two atoms have the same name
            if (atoms[i]->m_name == atoms[j]->m_name) {
                string error = "Each Atom in a Residue must have a unique name."
                               " These do not:\n";
                for(size_t k=0; k<atoms.size(); ++k) {
                    error += atoms[k]->str();}
                throw PANTZ_error (error);}}}
}
