/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the identification methods of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Search through a PDB file's information to identify the type of experiment
// that generated the data
void PROT::PDB::identify_type () {
    // Set the type of experiment to a default value
    m_type = "UNKNOWN";
    // This flag tells the program when to stop searching through the lines
    bool seen = false;
    // Loop through the lines
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with the appropriate phrase
        if (Text::startswith(m_lines[i], "EXPDTA")) {
            // Extract the phrase, going from the 11th character to the end of
            // the string
            string part = m_lines[i].substr(10);
            // If this is the first time the line has been seen
            if (!seen) {
                // Set the flag to true
                seen = true;
                // Set the type to the part
                m_type = part;}
            // If there are multiple lines of type information, concatenate
            // them appropriately
            else {
                // stick a semi-colon at the end of the current type
                // information if there is not one already
                if (m_type[m_type.size()-1] != ';') {m_type += ";";}
                m_type += " " + part;}}
        // If the experiment data type has already been stored, be done
        else if (seen) {break;}}
}

// Search through the PDB file's information for the resolution information
// for the experiment
void PROT::PDB::identify_resolution () {
    // Set the value to one that indicates no resolution information is
    // available
    m_resolution = -1.0;
    // Loop through the file's lines
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with the appropriate phrase
        if (Text::startswith(m_lines[i], "REMARK   2")) {
            // Split the line into pieces
            vector<string> parts = Text::split(m_lines[i]);
            // If there are the proper number of them
            if (parts.size() == 5) {
                // If the last part indicates a measurement
                if ((parts[4] == "ANGSTROMS.") || (parts[4] == "ANGSTROMS")) {
                    stringstream c1; c1 << parts[3]; c1 >> m_resolution;}
                // Stop the search
                break;}}}
}

// Use the REMARK 465 lines to identify missing residues
void PROT::PDB::identify_missing_residues (vector<vector<Residue> >& missing) {
    // Set up the missing vector to be accessible by character indexing
    missing.clear(); missing.resize(128);
    // This flag indicates when the relevant lines have been seen
    bool seen = false;
    // Store the appropriate lines in this vector
    vector<string> remarks;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with REMARK 465. I considered checking to see if
        // it was a REMARK line and then extracting the number, but this
        // SHOULD be more efficient in most cases.
        if (Text::startswith(m_lines[i], "REMARK 465")) {
            // Indicate that these lines have been seen
            seen = true;
            // If the line is long enough, store it with the header removed
            if (m_lines[i].size() > 10) {
                remarks.push_back(m_lines[i].substr(10));}}
        // Stop searching for lines once the REMARK 465 lines have been
        // collected or an ATOM line is seen
        else if (seen) {break;}
        else if (Text::startswith(m_lines[i], "ATOM")) {break;}}
    // If there were no REMARK 465 lines, be done
    if (remarks.size() == 0) {return;}
    // The standards for the REMARK 465 and 470 lines changed over time.
    // However, there is always an initial set of remarks and then a header
    // line that properly positions the contents. Find that header line
    size_t H = remarks.size();
    for(size_t i=0; i<remarks.size(); ++i) {
        if (Text::endswith(remarks[i], "RES C SSSEQI")) {
            H = i; break;}}
    // If that header line was not found, throw an error
    if (H == remarks.size()) {
        string text = "Failure to find the start of the Missing Residues.\n";
        throw PANTZ_error (text);}
    // Since the ending string is always the same, use math to determine the
    // appropriate columns for the information
    size_t n = remarks[H].size();
    // The column with the insertion character is n - 1
    size_t I = n-1;
    // The residue number can start 5 positions before that
    size_t N = I-5;
    // The chain is two positions before that
    size_t C = N-2;
    // And the name of the residue starts 4 positions before that
    size_t R = C-4;
    // Loop through the subsequent lines
    if (H == remarks.size() - 1) {return;}
    for (size_t i=H+1; i<remarks.size(); ++i) {
        // Get the name of the residues
        string resName = remarks[i].substr(R, 3); Text::strip(resName);
        // Get the name of the protein
        char L = remarks[i][C];
        // Get the number of the residue
        string numStr = remarks[i].substr(N, 5); Text::strip(numStr);
        long resNum = 0;
        stringstream c1; c1 << numStr; c1 >> resNum;
        // If there are enough characters, get the insertion code
        char insert = ' ';
        if (remarks[i].size() > I) {insert = remarks[i][I];}
        // Create a Residue using this information
        int j = (int) L;
        missing[j].push_back(Residue(resName, L));
        missing[j][missing[j].size()-1].m_number = resNum;
        missing[j][missing[j].size()-1].m_insertion = insert;
        missing[j][missing[j].size()-1].m_present = false;}
    // End the function
}

// The identify missing Atoms function is run after the Proteins have been
// assembled. It uses the REMARK 470 lines to store how many atoms are missing
// from various residues
void PROT::PDB::identify_missing_atoms () {
    // This flag indicates whether or not REMARK 470 lines have been seen
    bool seen = false;
    // Store the lines that list the missing Atoms here
    vector<string> remarks;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If this line is a REMARK 470 line
        if (Text::startswith(m_lines[i], "REMARK 470")) {
            // Indicate that the lines have been seen
            seen = true;
            // If the line is long enough, store it with the REMARK 470 part
            // removed
            if (m_lines[i].size() > 10) {
                remarks.push_back(m_lines[i].substr(10));}}
        // If the lines have been seen and this wasn't one, stop searching for
        // them
        else if (seen) {break;}
        // Or if an ATOM line is reached
        else if (Text::startswith(m_lines[i], "ATOM")) {break;}}
    // If there are no such remarks, end this function
    if (remarks.size() == 0) {return;}
    // The format for this section changed around 2007. The following code
    // should work correctly regardless of the format, but it takes a little
    // bit of work to correctly identify where the information is located.
    // First, find the line that ends the header information and describes
    // where the subsequent values are stored
    size_t H = remarks.size();
    for(size_t i=0; i<remarks.size(); ++i) {
        // In all cases, the line should end with CSSEQI  ATOMS
        if (Text::endswith(remarks[i], "CSSEQI  ATOMS")) {
            H = i; break;}}
    // If an appropriate line was not identified, throw an error
    if (H >= remarks.size() - 1) {
        string error = "Failure to identify missing atoms in REMARK 470 lines.\n";
        throw PANTZ_error (error);}
    // Get the length of that line
    size_t n = remarks[H].size();
    // The listing of the Atoms starts 4 positions before that
    size_t A = n - 4;
    // The insertion character is 4 positions before that
    size_t I = A - 4;
    // The residue number starts 4 positions before that
    size_t N = I - 4;
    // The chain starts 1 position before that
    size_t C = N - 1;
    // And the Residue starts either 4 or 5 positions before that
    size_t R = C - 5;
    if (remarks[H][C-4] == 'R') {R = C - 4;}
    // Loop through the lines listing missing Atoms
    for(size_t i=H+1; i<remarks.size(); ++i) {
        // Identify the residue's name
        string resName = remarks[i].substr(R, 3); Text::strip(resName);
        // Identify the Protein's name
        char L = remarks[i][C];
        // Extract the Residue's number's string
        string numStr = remarks[i].substr(N, 4); Text::strip(numStr);
        long resNum = 0;
        stringstream c1; c1 << numStr; c1 >> resNum;
        // Get the insertion character
        char insert = remarks[i][I];
        // Get the string of Atoms
        string atomStr = remarks[i].substr(A);
        // Split it into pieces
        vector<string> pieces = Text::split(atomStr);
        // A flag to indicate that the information was used
        bool flag = false;
        // Search the proteins to find the residue
        for(size_t j=0; j<m_proteins.size(); ++j) {
            // If this isn't the right protein, continue the search
            if (m_proteins[j].m_name != L) {continue;}
            // Check the protein's residues
            for(size_t k=0; k<m_proteins[j].m_count; ++k) {
                // Get a pointer to the Residue
                Residue * res = &(m_proteins[j].m_residues[k]);
                // If the Residue's information matches
                if ((res->m_name == resName) &&
                    ((res->m_number == resNum) &&
                     (res->m_insertion == insert))) {
                    // Indicate that the information is being used
                    flag = true;
                    // Add the number of missing atoms to the residue's
                    // information
                    res->m_missing_atoms += pieces.size();
                    break;}}
            // If the for loop reaches this point, it means it checked the
            // protein with the matching name. Break out of the search
            break;}
        // If the information was not used, throw an error
        if (!flag) {
            string error = "Failure to identify the Residue missing these "
                           "Atom:\n" + remarks[i] + "\n";
            throw PANTZ_error(error);}}
    // End the function
}
