/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to select a set of Atoms from the Residue. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Select a set of Atoms from the Residue
void PROT::Residue::select_atoms (vector<Atom *>& chosen, 
                                  const string how = "all") {
    // Make a copy of the provided selection string that's not a const
    string assess = how;
    // Make sure it is all lower case letters
    Text::lower(assess);
    // Convert the string to an integer method counter
    int method = 0;
    if (assess == "all") {method = 1;}
    else if (assess == "heavy") {method = 2;}
    else if (assess == "sidechain all") {method = 3;}
    else if (assess == "sidechain heavy") {method = 4;}
    else if (assess == "backbone all") {method = 5;}
    else if (assess == "backbone heavy") {method = 6;}
    else if (assess == "backbone main") {method = 7;}
    else if (assess == "CA") {method = 8;}
    else if (assess == "rotamer") {method = 9;}
    else {
        string error = how + " is not a recognized atom selection criteria.\n";
        throw PANTZ_error (error);}
    // if the Residue is empty, be done
    if (m_count == 0) {return;}
    // Implement the first two methods, which apply to any residue
    if (method == 1) {
        for(size_t i=0; i<m_count; ++i) {chosen.push_back(&(m_atoms[i]));}
        return;}
    if (method == 2) {
        for(size_t i=0; i<m_count; ++i) {
            if (!m_atoms[i].is_hydrogen()) {chosen.push_back(&(m_atoms[i]));}}
        return;}
    // The remaining methods only apply to amino acids, so error check that
    if (!is_amino_acid()) {
        string error = m_name + " is not an amino acid, so " + how + " is not "
                     "a valid atom selection criteria.\n";
        throw PANTZ_error (error);}
    // Side chain methods
    if (method == 3) {
        for(size_t i=0; i<m_count; ++i) {
            if (!m_atoms[i].is_backbone_atom()) {
                chosen.push_back(&(m_atoms[i]));}}
        return;}
    if (method == 4) {
        for (size_t i=0; i<m_count; ++i) {
            if ((!m_atoms[i].is_backbone_atom()) && (!m_atoms[i].is_hydrogen())) {
                chosen.push_back(&(m_atoms[i]));}}
        return;}
    // Backbone methods based on all atoms or only heavy atoms
    if (method == 5) {
        for (size_t i=0; i<m_count; ++i) {
            if (m_atoms[i].is_backbone_atom()) {
                chosen.push_back(&(m_atoms[i]));}}
        return;}
    if (method == 6) {
        for (size_t i=0; i<m_count; ++i) {
            if ((m_atoms[i].is_backbone_atom()) && (!m_atoms[i].is_hydrogen())) {
                chosen.push_back(&(m_atoms[i]));}}
        return;}
    // Access specific atom lists. These methods require all of the atoms to be
    // present.
    try {
        if (method == 7) {
            chosen.push_back(get_atom("N"));
            chosen.push_back(get_atom("CA"));
            chosen.push_back(get_atom("C"));
            return;}
        if (method == 8) {
            chosen.push_back(get_atom("CA"));
            return;}
        if (method == 9) {
            chosen.push_back(get_atom("CA"));
            chosen.push_back(get_atom("N"));
            if (m_name == "GLY") {chosen.push_back(get_atom("HA1"));}
            else {chosen.push_back(get_atom("CB"));}}}
    catch (PANTZ_error& e) {
        string error = "This error occurred when trying to select Atoms from a "
                       "Residue.\n";
        throw PANTZ_error (e, error);}
}

// Another version of the function that works with AtomPtr objects instead of
// pointers to Atoms
void PROT::Residue::select_atoms (vector<AtomPtr>& atoms, 
                                  const string how = "all") {
    // Make a vector of atom pointers
    vector<Atom *> chosen;
    // Select the atoms
    select_atoms(chosen, how);
    // Store the chosen pointers in the atoms vector
    for(size_t i=0; i<chosen.size(); ++i) {
        atoms.push_back(AtomPtr(chosen[i]));}
}

// Equivalent functions in the ResiduePtr class
void PROT::ResiduePtr::select_atoms (vector<Atom *>& atoms,
                                     const string how = "all") {
    check ();
    m_ptr->select_atoms(atoms, how);
}

void PROT::ResiduePtr::select_atoms (vector<AtomPtr>& atoms,
                                     const string how = "all") {
    check();
    m_ptr->select_atoms(atoms, how);
}
