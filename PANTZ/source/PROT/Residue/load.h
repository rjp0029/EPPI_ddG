/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * load method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The function that loads a Residue from a vector of Atom pointers
void PROT::Residue::load (const vector<PROT::Atom *>& atoms, 
                          const bool sidechain_only = false,
                          const bool complete_load = false) {
    // First, check that the Atoms are useable
    check_atoms (atoms);
    // If this is a complete load of the Residue's information, extract the
    // meaningful data from the atoms
    if (complete_load) {
        m_name = atoms[0]->m_residue;
        m_number = atoms[0]->m_residue_number;
        m_insertion = atoms[0]->m_insertion;
        m_protein = atoms[0]->m_protein;
        // Set the internal number to 0 - it will have to be externally
        // specified
        m_internal = 0;
        // Because the Residue is being loaded from Atoms, it is definitely
        // present
        m_present = true;}
    // If a complete load is not being done and a side chain load is not being
    // done, make sure the provided atoms have the proper residue name
    if ((!complete_load) && (!sidechain_only)) {
        if (atoms[0]->m_residue != m_name) {
            string error = "It is not permitted to load a "
                         + atoms[0]->m_residue + " in place of a "
                         + m_name + ".\n";
            throw PANTZ_error (error);}}
    // The relevant Atoms have to be collected. Store them here
    vector<Atom> use;
    // If a side-chain only load is being done
    if (sidechain_only) {
        // If the residue is not an amino acid, throw an error
        if (!CHECK::is_amino_acid(m_name)) {
            string error = m_name + " is not an amino acid, so a "
                         "sidechain-only Residue load is not allowed.\n";
            throw PANTZ_error (error);}
        // Update the name of the residue to match that from the provided
        // atoms
        m_name = atoms[0]->m_residue;
        // Loop through the current atoms to identify the backbone atoms
        for(size_t i=0; i<m_count; ++i) {
            if (CHECK::is_backbone_atom(m_atoms[i].m_name)) {
                // Make sure proper conventions are followed for glycine and
                // proline
                if ((m_name == "GLY") && (m_atoms[i].m_name == "HA")) {
                    m_atoms[i].m_name = "HA2";}
                else if ((m_name != "GLY") && (m_atoms[i].m_name == "HA2")) {
                    m_atoms[i].m_name = "HA";}
                else if (m_name == "PRO") {
                    if (m_atoms[i].m_name == "HN") {continue;}
                    else if (m_atoms[i].m_name == "HT3") {continue;}
                    else if (m_atoms[i].m_name == "HT1") {
                        m_atoms[i].m_name = "HN1";}
                    else if (m_atoms[i].m_name == "HT2") {
                        m_atoms[i].m_name = "HN2";}}
                else if ((m_atoms[i].m_name == "HN1") && (m_name != "PRO")) {
                    m_atoms[i].m_name = "HT1";}
                else if ((m_atoms[i].m_name == "HN2") && (m_name != "PRO")) {
                    m_atoms[i].m_name = "HT2";}
                use.push_back(m_atoms[i]);}}
        // Loop through the provided atoms to identify the non-backbone atoms
        for(size_t i=0; i<atoms.size(); ++i) {
            if (!CHECK::is_backbone_atom(atoms[i]->m_name)) {
                use.push_back(*(atoms[i]));}}}
    else {
        for(size_t i=0; i<atoms.size(); ++i) {
            use.push_back(*(atoms[i]));}}
            // Otherwise, use all of the atoms
    // Delete existing atoms from the residue
    clean_up();
    // Allocate space for the new set of atoms
    m_count = use.size();
    m_atoms = new Atom [m_count];
    // Store the chosen atoms
    for(size_t i=0; i<m_count; ++i) {m_atoms[i] = use[i];}
    // Make sure that every atom has the proper labelling information
    for(size_t i=0; i<m_count; ++i) {
        m_atoms[i].m_alt = ' ';
        m_atoms[i].m_residue = m_name;
        m_atoms[i].m_residue_number = m_number;
        m_atoms[i].m_insertion = m_insertion;
        m_atoms[i].m_protein = m_protein;}
    // Ensuring consistent Atom numbering is done outside of this function. End
    // the function.
}

// Load a Residue from a vector of AtomPtrs instead of a vector of Atom Pointers
void PROT::Residue::load (vector<AtomPtr>& atoms, 
                          const bool sidechain_only = false,
                          const bool complete_load = false) {
    // Make a vector of just the pointers
    vector<Atom *> ptrs; ptrs.reserve(atoms.size());
    for(size_t i=0; i<atoms.size(); ++i) {ptrs.push_back(atoms[i].pointer());}
    // Load them
    load(ptrs, sidechain_only, complete_load);
}
