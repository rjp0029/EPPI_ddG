/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the residue initialization and integration methods of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Use the SEQRES information to set up initial vectors of Residues
void PROT::PDB::initialize_Residues (vector<vector<Residue> >& residues) {
    // Set up the residues so that the vector can directly access contents via
    // character indexing
    residues.clear(); residues.resize(128);
    // This flag indicates when the necessary lines have been observed
    bool seen = false;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line startswith SEQRES
        if (Text::startswith(m_lines[i], "SEQRES")) {
            seen = true;
            // The chain must be in character 12
            char L = m_lines[i][11];
            // Extract the substring containing the residues (positions 20 -
            // end of string)
            string res_part = m_lines[i].substr(19);
            // strip whitespace
            Text::strip(res_part);
            // Split into pieces
            vector<string> parts = Text::split(res_part);
            // Create and store each of them as a Residue
            for(size_t j=0; j<parts.size(); ++j) {
                residues[(int) L].push_back(Residue(parts[j], L));}}
        // If the sequence residue information has been found, be done
        else if (seen) {break;}}
}

// Integrate Missing Residues into the list of Residues that make up a
// Protein. Then store that Protein in the PDB file's internal vector
void PROT::PDB::integrate_missing_residues (Protein& prot,
                                            vector<Residue>& missing) {
    // If there are no missing Residues, just store the protein as it is
    if (missing.size() == 0) {m_proteins.push_back(prot); return;}
    // Since there are missing Residues, extract the Residues from the protein
    // and store it as a vector
    vector<Residue> residues; residues.reserve(prot.m_count + missing.size());
    for(size_t i=0; i<prot.m_count; ++i) {
        residues.push_back(prot.m_residues[i]);}
    // Loop through each missing residue
    for(size_t i=0; i<missing.size(); ++i) {
        // A pointer to the missing residue
        Residue * m = &(missing[i]);
        // The vector of residues needs to be searched for the best place to
        // insert this missing residue. Store the best insertion score and the
        // index of the residue this one should be placed before here
        size_t best_score = residues.size() * 100000;
        size_t best_index = residues.size();
        // A pointer to a "Known" residue, starting with the first
        Residue * k = &(residues[0]);
        // Determine if the missing residue goes before the first residue in
        // the protein
        if ((m->m_number < k->m_number) || 
           ((m->m_number == k->m_number) && (m->m_insertion < k->m_insertion))){
            // If it does, the score is guaranteed to be better than the best
            // score, so just save that information
            best_score = (k->m_number - m->m_number) + 1;
            best_index = 0;}
        // Check to see if the missing residue inserts well into other points
        // in the vector of Residues
        for(size_t j=1; j<residues.size(); ++j) {
            // The known residue is this one
            k = &(residues[j]);
            // And a previous residue is the one before it
            Residue * p = &(residues[j-1]);
            // The missing residue must come after the previous residue, so if
            // it doesn't continue the search
            if ((m->m_number < p->m_number) ||
               ((m->m_number == p->m_number) && 
                (m->m_insertion < p->m_insertion))) {continue;}
            // The missing residue must come before the current residue
            if ((k->m_number < m->m_number) ||
               ((k->m_number == m->m_number) && 
                (k->m_insertion < m->m_insertion))) {continue;}
            // Since the missing residue comes after the previous residue and
            // before the current residue, calculate a score for inserting it
            // here
            // This calculation is weird because I've been able to identify at
            // least one case where units of thousands are used to indicate
            // insertions instead of an insertion character
            size_t score;
            if (k->m_number >= p->m_number) {score = k->m_number - p->m_number;}
            else {score = p->m_number - k->m_number;}
            // If this is the best score found so far
            if (score < best_score) {
                best_score = score;
                best_index = j;}}
        // Calculate a score for inserting the residue at the end of the
        // vector
        k = &(residues[residues.size()-1]);
        if ((m->m_number > k->m_number) ||
           ((m->m_number == k->m_number) && (m->m_insertion > k->m_insertion))){
            size_t score = (m->m_number - k->m_number) + 1;
            if (score < best_score) {
                best_score = score;
                best_index = residues.size();}}
        // Insert the missing residue into the best position identified for
        // it. If that is the end of the vector, just use push_back to store
        // it
        if (best_index == residues.size()) {residues.push_back(missing[i]);}
        // Otherwise, insert the residue at the proper position. To do that, I
        // need an iterator to the start of the vector
        else {
            vector<Residue>::iterator it = residues.begin();
            // Insert the residue
            residues.insert(it + best_index, missing[i]);}}
    // Store the residues as a Protein in the PDB file
    m_proteins.push_back(Protein(residues));
    // End the function
}
