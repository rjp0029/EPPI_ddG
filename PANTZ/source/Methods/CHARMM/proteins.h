/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains functions that deal with outputting and reading in
 * Proteins for CHARMM calculations */

// This file should only be included by CHARMM.h
#ifndef CHARMM_Loading_Status
#error CHARMM methods must be included by CHARMM.h
#endif

// Prepare proteins for CHARMM calculations and tell CHARMM to load them
void CHARMM::proteins_for_charmm (string& script,
                                  vector<PROT::Protein>& prots) {
    // Atom and residue numbering information
    long rn = 1;
    long an = 1;
    // Do this for each protein
    for(size_t i=0; i<prots.size(); ++i) {
        // Renumber the protein's residues and atoms sequentially
        rn = prots[i].renumber_residues(rn);
        an = prots[i].renumber_atoms(an);
        // change the names of histidine residues to HSD
        prots[i].for_charmm_histidine_fix ();
        // Get the name of a file for this protein
        string fileName = make_protein_name (prots[i], true);
        // Write the protein to that file
        ofstream output; output.open(fileName.c_str());
        if (!output.is_open()) {
            string error = "Failed to open " + fileName + " for charmm\n";
            throw PANTZ_error (error);}
        // Write hte protein to the file
        output << prots[i].str(true) << "END\n";
        output.close();
        // Update the script with how to load this protein
        script += "! Load Protein ";
        script += prots[i].name();
        script += "\nopen read unit 10 form name " + fileName + "\n"
                  "read sequ pdb offi unit 10\n"
                  "close unit 10\n"
                  "gene pr";
        char L = prots[i].name();
        string name = "";
        name += L;
        Text::lower(name);
        script += name + " setup\n"
                  "open read unit 10 form name " + fileName + "\n"
                  "read coor pdb unit 10\n"
                  "close unit 10\n\n";}
}

// Load proteins after charmm energy calculations
void CHARMM::proteins_from_charmm (vector<PROT::Protein>& prots) {
    // Do this for each protein
    for(size_t i=0; i<prots.size(); ++i) {
        // Get the name of the file
        string fileName = make_protein_name (prots[i], false);
        // Load the protein from that file
        prots[i].load(fileName);
        // Fix histidines
        prots[i].from_charmm_histidine_fix ();}
}
