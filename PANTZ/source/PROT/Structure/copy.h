/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Structure.h
 * header file, and includes pre-processor directives for that behavior. The
 * copy methods of the Structure class are defined here. */

// Make sure the Structure class is currently being loaded
#ifndef Structure_Loading_Status
#error Structure methods must be included by the Structure.h header file
#endif

// Copy contents from another Structure
void PROT::Structure::copy (const Structure * other) {
    // Copy the protein pointers
    m_proteins.clear(); 
    if (other->m_proteins.size() > 0) {
        m_proteins.reserve(other->m_proteins.size());
        for(size_t i=0; i<other->m_proteins.size(); ++i) {
            m_proteins.push_back(other->m_proteins[i]);}}
    // Copy the names
    m_names.clear();
    if (other->m_names.size() > 0) {
        m_names.reserve(other->m_names.size());
        for(size_t i=0; i<other->m_names.size(); ++i) {
            m_names.push_back(other->m_names[i]);}}
}

// Copy the contents of another Structure, but create pointers to a different
// list of proteins. This is used in copying PDB files
void PROT::Structure::copy (const Structure * other, vector<Protein>& prots) {
    // Copy the names of the other structure - this copies the names and sets up
    // the proteins vector to be the right size
    copy(other);
    // If there are proteins
    if(m_proteins.size() > 0) {
        // Loop through them
        for(size_t i=0; i<m_proteins.size(); ++i) {
            // Get the name of the current protein
            char L = m_proteins[i]->name();
            // Use a boolean flag to make sure that protein is found
            bool flag = false;
            // Find that protein in the provided vector
            for(size_t j=0; j<prots.size(); ++j) {
                if (prots[j].name() == L) {
                    flag = true;
                    m_proteins[i] = &(prots[j]);
                    break;}}
            // If the protein wasn't found, throw an error
            if (!flag) {
                string error = "Failure to find Protein ";
                error += L;
                error += " by the Structure copy function.\n";
                throw PANTZ_error (error);}}}
}
