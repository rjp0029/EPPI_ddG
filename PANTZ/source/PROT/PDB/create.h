/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the Protein construction and structure creation methods of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Construct the Proteins from the information in the PDB file
void PROT::PDB::construct_Proteins () {
    // Collect the information on the SEQRES data
    vector<vector<Residue> > seqres;
    initialize_Residues(seqres);
    // Collect the information on the missing Residues
    vector<vector<Residue> > missing;
    identify_missing_residues(missing);
    // Collect the Atoms that make up the Proteins
    vector<vector<Atom> > atoms;
    collect_Atoms (atoms);
    // Generate proteins with those Atoms
    vector<Protein> prots;
    for(size_t i=0; i<128; ++i) {
        if (atoms[i].size() > 0) {prots.push_back(Protein(atoms[i]));}}
    // If there are no proteins, raise an error
    if (prots.size() == 0) {
        string error = "No Proteins were identified.\n";
        throw PANTZ_error (error);}
    // Loop through those proteins
    for(size_t P=0; P<prots.size(); ++P) {
        // Get the character index of this protein
        int i = (int) prots[P].m_name;
        // Integrate the missing residues into that protein and store it in
        // the PDB file
        integrate_missing_residues (prots[P], missing[i]);
        // Validate the SEQRES information for this protein if there is SEQRES
        // data
        if (seqres[i].size() > 0) {
            validate_seqres (&(m_proteins[m_proteins.size()-1]), seqres[i]);}}
    // Identify atoms that may be missing in those proteins
    identify_missing_atoms ();
}

// After the Proteins have been constructed, they can be grouped into
// Structures
void PROT::PDB::create_Structures () {
    // Identify the compound line statements and separate them based on the
    // molecule ID
    vector<string> compound;
    vector<vector<string> > compounds;
    bool seen = false;
    // Loop through the lines of the PDB file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If this is the right sort of line
        if (Text::startswith(m_lines[i], "COMPND")) {
            seen = true;
            // Extract the useful part of the line
            string phrase = m_lines[i].substr(10);
            Text::strip(phrase);
            // If this is a Molecule ID line and there is existing compound
            // information, store it
            if (Text::startswith(phrase, "MOL_ID:")) {
                if (compound.size() > 0) {
                    compounds.push_back(compound);
                    compound.clear();}}
            // Otherwise, store this line in the compound
            else {compound.push_back(phrase);}}
        else if (seen) {break;}}
    if (compound.size() > 0) {compounds.push_back(compound);}
    // If no compound lines were seen, end this function
    if (!seen) {return;}
    // Create a vector of boolean values indicating which Proteins have been
    // assigned to a Structure
    vector<bool> used; used.reserve(m_proteins.size());
    for(size_t i=0; i<m_proteins.size(); ++i) {used[i] = false;}
    // Loop through the compound information
    for(size_t C=0; C<compounds.size(); ++C) {
        // Store a structure in the structures vector
        m_structures.push_back(Structure());
        size_t n = m_structures.size() - 1;
        // Loop through the relevant lines
        for(size_t i=0; i<compounds[C].size(); ++i) {
            // I only care about some of the lines. For instance, if the line
            // lists one or more names of the structure
            if ((Text::startswith(compounds[C][i], "MOLECULE:")) ||
                (Text::startswith(compounds[C][i], "SYNONYM:"))) {
                // Extract the meaningful part of the line and strip it of
                // whitespace
                string text = compounds[C][i].substr(9);
                Text::strip(text);
                // Split the text on commas, which separate the names
                vector<string> parts = Text::split(text, ',');
                // Loop through the parts
                for(size_t j=0; j<parts.size(); ++j) {
                    // Strip the part of whitespace
                    Text::strip(parts[j]);
                    // If the last character is a semi-colon, delete it
                    if (parts[j][parts[j].size()-1] == ';') {
                        parts[j].erase(parts[j].size()-1, 1);}
                    // Store the name
                    m_structures[n].m_names.push_back(parts[j]);}}
            else if (Text::startswith(compounds[C][i], "CHAIN:")) {
                // Strip off the leading portion of the string and whitespace
                string text = compounds[C][i].substr(6);
                Text::strip(text);
                // Split it on commas
                vector<string>parts = Text::split(text, ',');
                // Loop through that list
                for(size_t j=0; j<parts.size(); ++j) {
                    // Strip whitespace from the string
                    Text::strip(parts[j]);
                    // Get the first character
                    char L = parts[j][0];
                    // Use this flag to indicate whether or not that protein
                    // is found
                    bool found = false;
                    for(size_t k=0; k<m_proteins.size(); ++k) {
                        if (m_proteins[k].m_name == L) {
                            found = true;
                            // If this protein has already been used, raise an
                            // error
                            if (used[k]) {
                                string error = "Chain ";
                                error += L;
                                error += "is used in more than one Molecule.\n";
                                throw PANTZ_error (error);}
                            // Indicate that the protein is being used
                            used[k] = true;
                            // Store the pointer in the structure
                            m_structures[n].m_proteins.push_back(&(m_proteins[k]));
                            break;}}
                    // If the protein wasn't identified, raise an error
                    if (!found) {
                        string error = "Chain ";
                        error += L;
                        error += " is listed in the COMPND lines, but no Atoms "
                                 "were identified for it.\n";
                        throw PANTZ_error (error);}}}}}
    // If there are any proteins that have not been used, store them in an
    // unamed structure
    bool flag = false;
    size_t n = m_structures.size();
    for(size_t i=0; i<m_proteins.size(); ++i) {
        if (!used[i]) {
            if (!flag) {
                flag = true;
                m_structures.push_back(Structure());}
            m_structures[n].m_proteins.push_back(&(m_proteins[i]));}}
    // End the function
}
