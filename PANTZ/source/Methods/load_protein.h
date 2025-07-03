/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the protein loading method. */

// This file is supposed to be loaded by Methods.h
#ifndef Methods_Loading_Status
#error General Methods must be included by Methods.h
#endif

// Implement the load protein function
void METHODS::load_protein (vector<string> * command,
                            vector<PROT::Protein>& proteins, 
                            vector<PROT::PDB>& pdbs,
                            ofstream& output) {
    // Validate that the command has 5 entries
    if (command->size() != 5) {
        string error = "Algorithm error: METHODS::load_protein called with "
                       "invalid command\n";
        throw PANTZ_error (error);}
    // Get information from the command. The first entry is a classification
    // that this function actually doesn't have to worry about - it's used where
    // it is called from to change the provided proteins vector
    string path = (*command)[1];
    string fileName = (*command)[2];
    char L1 = (*command)[3][0];
    char L2 = (*command)[4][0];
    // Confirm that no existing protein is named L2
    for(size_t i=0; i<proteins.size(); ++i) {
        if (proteins[i].name() == L2) {
            string error = "Chain ";
            error += L1;
            error += " from " + path + fileName + " could not be loaded as ";
            error += L2;
            error += " because a protein with that name is already loaded.\n";
            throw PANTZ_error (error);}}
    // Whether or not the PDB files that were already loaded have this protein
    // in them
    bool known = false;
    size_t index = 0;
    // Go through the existing PDB files
    for(size_t i=0; i<pdbs.size(); ++i) {
        // If the PDB file's folder and name match
        if ((pdbs[i].folder() == path) && (pdbs[i].name() == fileName)) {
            known = true; index = i; break;}}
    // If the file is not already known, load it
    if (!known) {
        index = pdbs.size();
        pdbs.push_back(PROT::PDB(fileName, path));}
    // Find the specified protein
    known = false;
    for(size_t i=0; i<pdbs[index].proteins(); ++i) {
        // Get the pointer
        PROT::Protein * prot = pdbs[index].protein(i);
        // If it is the right protein
        if (prot->name() == L1) {
            // Mark it as found
            known = true;
            // Make a vector of Residues from the protein where the residues are
            // amino acids that are complete
            vector<PROT::Residue> residues;
            // Go through the protein's residues
            for(size_t j=0; j<prot->size(); ++j) {
                // Get the Residue's pointer
                PROT::ResiduePtr res;
                res = (*prot)[j];
                // if it is not an amino acid, skip it
                if (!res.is_amino_acid()) {continue;}
                // If it's score is > 0, use it
                if (res.score () > 0) {
                    // Store a duplicate of the residue
                    residues.push_back(res.duplicate());}}
            // If no residues were found
            if (residues.size() == 0) {
                string error = "Chain ";
                error += L1;
                error += " from " + path + fileName + " did not contain "
                         "any amino acids.\n";
                throw PANTZ_error(error);}
            // Store it
            proteins.push_back(PROT::Protein(residues));
            // Rename it
            proteins[proteins.size()-1].set_name(L2);
            // Stop the search
            break;}}
    // If it wasn't found, throw an error
    if (!known) {
        string error = "Chain ";
        error += L1;
        error += " was not found in " + path + fileName + "\n";
        throw PANTZ_error (error);}
    // if a summary file exists, add a message to it
    if (output.is_open()) {
        string text = "Protein ";
        text += L2;
        text += " loaded on " + METHODS::time_stamp() + "\n";
        output << text << endl;}
    // End the function
}
