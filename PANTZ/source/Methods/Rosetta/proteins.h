/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements protein I/O functions for Rosetta calculations */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// Output proteins for use in Rosetta calculations
void Rosetta::proteins_for_rosetta (vector<PROT::Protein>& proteins,
              const string& fileLabel, const bool internal = true) {
    // Setup atom and residue number values
    long resNum = 1;
    long atomNum = 1;
    if (!internal) {resNum = 0;}
    // Open the appropriate file
    string fileName = fileLabel + ".pdb";
    ofstream output; output.open(fileName.c_str());
    // If it failed to open, throw an error
    if (!output.is_open()) {
        string error = "Failed to open " + fileName + " for Rosetta "
                       "calculations.\n";
        throw PANTZ_error (error);}
    // Write the content of each protein to that file
    for(size_t i=0; i<proteins.size(); ++i) {
        // The rosetta_str method of the protein class (really of the residue
        // class) takes care of the numbering requirements, stripping out
        // hydrogens, renaming atoms (CD of ILE to CD1), and renumbering the
        // residues and atoms as appropriate
        output << proteins[i].rosetta_str(resNum, atomNum);}
    // Put an 'end' at the end of the file
    output << "END\n";
    output.close();
}

// Read in proteins from Rosetta calculations
void Rosetta::proteins_from_rosetta (vector<PROT::Protein>& proteins,
                                     const string& fileLabel) {
    // Create a vector of vector of Atoms
    vector<vector<PROT::Atom> > allAtoms;
    // And a vector that is those in a particular protein
    vector<PROT::Atom> atoms;
    // The character of the protein of the most recent atom assessed
    char L = 0;
    // The contents should be in this file
    string fileName = fileLabel + "_0001.pdb";
    // REad them in
    ifstream input; input.open(fileName.c_str());
    if (!input.is_open()) {
        string error = "Failed to open " + fileName + " after Rosetta "
                       "calculations.\n";
        throw PANTZ_error (error);}
    // Get the first line
    string line; getline(input, line);
    // Go through the file
    while (!input.eof()) {
        // If the line starts with ATOM
        if (Text::startswith(line, "ATOM")) {
            // Make an Atom
            PROT::Atom atom (line);
            // If the last protein character isn't known, set it equal to this
            // one
            if (L == 0) {L = atom.protein();}
            // If this atom is part of a new protein, update the atom storage
            if (L != atom.protein()) {
                if (atoms.size() > 0) {
                    allAtoms.push_back(atoms);
                    atoms.clear();}
                L = atom.protein();}
            // Store the atom
            atoms.push_back(atom);}
        // Get the next line
        getline(input, line);}
    // If there are unstored atoms, store them
    if (atoms.size() > 0) {allAtoms.push_back(atoms); atoms.clear();}
    // Close the input file
    input.close();
    // If the number of vectors in all atoms doesn't match the number of
    // proteins, throw an error
    if (allAtoms.size() != proteins.size()) {
        stringstream c1; c1 << allAtoms.size();
        stringstream c2; c2 << proteins.size();
        string error = fileName + " contained " + c1.str() + " protein chains "
                       "after Rosetta calculations when " + c2.str() + " were "
                       "input.\n";
        throw PANTZ_error (error);}
    // Load each protein from the appropriate set of atoms
    for(size_t i=0; i<proteins.size(); ++i) {
        // Whether or not the protein has been found and loaded
        bool loaded = false;
        for(size_t j=0; j<allAtoms.size(); ++j) {
            if (allAtoms[j][0].protein() == proteins[i].name()) {
                // Indicate the protein was found
                loaded = true;
                // Load the protein from that vector of Atoms
                proteins[i].load(allAtoms[j]);
                // Update the atom names after Rosetta
                proteins[i].update_atoms_after_Rosetta();
                // stop the search
                break;}}
        // If the protein wasn't loaded, throw an error
        if (!loaded) {
            string error = "Protein ";
            error += proteins[i].name();
            error += " not found in " + fileName + "\n";
            throw PANTZ_error (error);}}
    // End the function
}
