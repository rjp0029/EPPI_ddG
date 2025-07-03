/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the check_residues method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Determine whether or not a vector of Residues is acceptable for use in a
// Protein
void PROT::Protein::check_residues (vector<PROT::Residue>& residues) const {
    // If there are no residues, raise an error
    if (residues.size() == 0) {
        string text = "A Protein cannot be loaded from an empty vector of "
                      "Residues.\n";
        throw PANTZ_error (text);}
    // If there is only a single residue, be done 
    else if (residues.size() == 1) {return;}
    // Since there are multiple residues, make sure they all have the same
    // protein character
    for(size_t i=1; i<residues.size(); ++i) {
        if (residues[0].protein() != residues[i].protein()) {
            string text = "A Protein cannot be loaded from a vector of "
                          "Residues that belong to different proteins.\n";
            throw PANTZ_error (text);}}
    // Check the residues to see which (if any) of them have only HETATMs
    vector<bool> hetatms; hetatms.reserve(residues.size());
    for(size_t i=0; i<residues.size(); ++i) {
        bool het = true;
        for(size_t j=0; j<residues[i].size(); ++j) {
            if (residues[i].m_atoms[j].type() != "HETATM") {
                het = false; break;}}
        hetatms.push_back(het);}
    // Ensure every Residue has a unique number / insertion code combination
    for(size_t i=0; i<residues.size() - 1; ++i) {
        // Skip this residue if it only contains hetero atoms
        if (hetatms[i]) {continue;}
        // Loop through all subsequent residues
        for (size_t j=i+1; j<residues.size(); ++j) {
            if (hetatms[j]) {continue;}
            // Check the residues' numbering data
            if ((residues[i].number() == residues[j].number()) &&
                (residues[i].insertion_code() == residues[j].insertion_code())) {
                string error = "Every Residue in a Protein must have a unique "
                               "number.\n";
                error += residues[i].str() + residues[j].str();
                throw PANTZ_error (error);}}}
}
