/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the SEQRES validation method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Validate that identified SEQRES information matches the specified Protein's
// information
void PROT::PDB::validate_seqres (Protein * prot, vector<Residue>& seqres) {
    // There must be at least as many residues in the protein as there are in
    // the sequence residue data
    if (prot->m_count < seqres.size()) {
        string error = "The SEQRES data contains more residues than are in the "
                       "ATOM and Missing Residue lines for chain ";
        error += prot->m_name;
        throw PANTZ_error (error + ".\n");}
    // Create a vector of boolean values indicating whether or not the a
    // protein residue should be compared to the SEQRES data
    vector<bool> use; use.resize(prot->m_count);
    // Go through every protein residue
    for(size_t i=0; i<prot->m_count; ++i) {
        // Initially set the value to false (don't use)
        use[i] = false;
        // As a short-cut for HETATM entries like water
        if ((i > 0) && 
            (prot->m_residues[i].m_name == prot->m_residues[i-1].m_name)) {
            use[i] = use[i-1]; continue;}
        // Search the SEQRES data
        for(size_t j=0; j<seqres.size(); ++j) {
            // If there is at least one SEQRES entry with the same Residue
            // name, use this residue in the comparison
            if (seqres[j].m_name == prot->m_residues[i].m_name) {
                use[i] = true; break;}}}
    // Count up the number of usable residues
    size_t count = 0;
    for(size_t i=0; i<use.size(); ++i) {if(use[i]) {++count;}}
    // This count MUST match the number of SEQRES entries
    if (count != seqres.size()) {
        stringstream c1; c1 << seqres.size();
        stringstream c2; c2 << count;
        string error = "Inconsistency between the SEQRES and Sequence "
                       "information for chain ";
        error += prot->m_name;
        error += ".\nThe SEQRES data lists " + c1.str() + " residues "
                 "but " + c2.str() + " residues have sequence-eligible "
                 "structure information.\n";
        throw PANTZ_error (error);}
    // An index for which residue to use in the Protein
    size_t PI = 0;
    // Loop through the SEQRES values
    for(size_t i=0; i<seqres.size(); ++i) {
        // Find the next protein residue to use
        while (!use[PI]) {++PI;}
        // If the residue information does not match, raise an error
        if (prot->m_residues[PI].m_name != seqres[i].m_name) {
            stringstream c1; c1 << i+1;
            stringstream c2; c2 << PI+1;
            stringstream c3; c3 << prot->m_residues[PI].m_number;
            string resNum = c3.str();
            resNum += prot->m_residues[PI].m_insertion;
            Text::strip(resNum);
            string error = "A Mismatch occurred between the SEQRES and "
                           "identified sequence information in chain ";
            error += prot->m_name;
            error += ":\nThis happend at residue " + c1.str() + " in the "
                     "SEQRES data and residue " + c2.str() + " (" + resNum
                   + ") in the sequence information.\n";
            throw PANTZ_error (error);}
        // Increment the protein index
        ++PI;}
}
