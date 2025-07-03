/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the load methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Load a Protein from a vector of residues
void PROT::Protein::load (vector<Residue>& residues) {
    // First, check that the provided residues are useable. This does several
    // things, including make sure there is at least one residue
    check_residues (residues);
    // If the protein is empty, do one thing
    if (m_residues == 0) {
        // Store the number of residues
        m_count = residues.size();
        // Allocate memory
        m_residues = new Residue [m_count];
        // Store the name of the protein
        m_name = residues[0].protein();
        // Store the residues and update their internal numbers as they are
        // stored
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i] = residues[i];
            m_residues[i].m_internal = i+1;}}
    // If the protein has previous content
    else {
        // if the number of residues doesn't match
        if (residues.size() != m_count) {
            string error = "A previously-loaded Protein can only be reloaded "
                           "from a vector with the same number of Residues as "
                           "the Protein already contains.\n";
            stringstream c1; c1 << m_count;
            stringstream c2; c2 << residues.size();
            error += "Protein ";
            error += m_name;
            error += " contains " + c1.str() + " residues. It is being "
                     "reloaded from a vector of " + c2.str() + " residues.\n";
            throw PANTZ_error (error);}
        // Loop through the residues
        for(size_t i=0; i<m_count; ++i) {
            // Make a vector of pointers to the Residue's atoms
            vector<Atom *> ptrs; ptrs.reserve(residues[i].size());
            residues[i].select_atoms(ptrs, "all");
            // Load them into the residue of the protein
            m_residues[i].load(ptrs);}}
    // set the terminus status of the first and last residues
    m_residues[0].m_N_terminus = true;
    m_residues[m_count-1].m_C_terminus = true;
    // make sure atoms are sequentially numbered
    long n = renumber_atoms (1);
}

// Load a Protein from a vector of Atoms
void PROT::Protein::load (vector<Atom>& atoms) {
    // Make sure there are atoms
    if (atoms.size() == 0) {
        string error = "A Protein cannot be loaded from an empty vector of "
                       "Atoms.\n";
        throw PANTZ_error (error);}
    // Store pointers to Atoms in the same residue in this vector
    vector<Atom *> res;
    // The labelling information to distinguish residues
    long n = atoms[0].residue_number();
    char l = atoms[0].insertion_code();
    // Store the atoms as residues
    vector<Residue> residues;
    // Loop through the Atoms
    for(size_t i=0; i<atoms.size(); ++i) {
        // If this is a different residue than the previous one
        if ((atoms[i].residue_number() != n) ||
            (atoms[i].insertion_code() != l)) {
            // Store the residue (if it is not empty, which it should never be)
            if (res.size() > 0) {residues.push_back(Residue(res));}
            // Clear that vector
            res.clear();
            // Update the labelling information
            n = atoms[i].residue_number();
            l = atoms[i].insertion_code();}
        // Only store the atom if it's alternative location character is blank
        // or A
        if ((atoms[i].alternative_location() == ' ') || 
            (atoms[i].alternative_location() == 'A')) {
            res.push_back(&(atoms[i]));}}
    // If the res isn't empty, store it
    if (res.size() > 0) {residues.push_back(res);}
    // Load the residues
    try {load(residues);}
    catch (PANTZ_error& e) {
        string error = "This error occured in the Atom-based load method of "
                       "the Protein class.\n";
        throw PANTZ_error (e, error);}
}

// Load a Protein from a vector of strings
void PROT::Protein::load (const vector<string>& contents) {
    // make sure there are values
    if (contents.size() == 0) {
        string error = "A Protein cannot be loaded from an empty vector of "
                       "strings.\n";
        throw PANTZ_error (error);}
    // Convert them to Atoms
    vector<Atom> atoms;
    for(size_t i=0; i<contents.size(); ++i) {
        // not everything will might be an Atom, so only store the Atom objects
        try{atoms.push_back(Atom(contents[i]));}
        catch(PANTZ_error) {continue;}}
    // Load them
    try {load(atoms);}
    catch (PANTZ_error& e) {
        string error = "This error occurred in the string-based load method of "
                       "the Protein class.\n";
        throw PANTZ_error (e, error);}
}

// Load a Protein from a file
void PROT::Protein::load (const string& fileName, const string path = "./") {
    // create the file name for loading
    string use = path + fileName;
    // Attempt to open that file
    ifstream input; input.open(use.c_str());
    if (!input.is_open()) {
        string error = "Failure to open this Protein file:\nName: "
                     + fileName + "\nLocation: " + path + "\n";
        throw PANTZ_error(error);}
    // Store the contents as string
    vector<string> contents;
    string line; getline(input, line);
    while (!input.eof()) {
        contents.push_back(line); getline(input, line);}
    input.close();
    // Load them
    try {load(contents);}
    catch (PANTZ_error& e) {
        string error = "This error occurred in the file-based load method of "
                       "the Protein class.\n";
        throw PANTZ_error (e, error);}
}

// Constructing a Protein using two strings
PROT::Protein::Protein (const string& fileName, const string path = "./") {
    initialize(); load(fileName, path);
}
