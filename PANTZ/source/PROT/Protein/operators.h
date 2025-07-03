/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the operator methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Provide access to the Protein's residues
PROT::Residue * PROT::Protein::operator() (const long num,
                                           const char insertion = ' ',
                                           const bool position = false) const {
    // Throw an error if the Protein is empty
    if (m_count == 0) {
        string text = "It is not possible to access a Residue in an empty "
                      "Protein.\n";
        throw PANTZ_error (text);}
    // If the residue is being referenced by position
    if (position) {
        // If the number is in the proper range
        if ((num >= 0) && (num < m_count)) {
            return &(m_residues[num]);}
        else {
            stringstream c1; c1 << num;
            stringstream c2; c2 << m_count;
            string error = c1.str() + " is not a valid Residue index in a "
                           "Protein with " + c2.str() + " Residues.\n";
            throw PANTZ_error (error);}}
    // Otherwise, check the residues for number / insertion code combinations
    // that match the provided values
    else {
        for(size_t i=0; i<m_count; ++i) {
            if ((m_residues[i].m_number == num) &&
                (m_residues[i].m_insertion == insertion)) {
                return &(m_residues[i]);}}
        stringstream c1; c1 << num;
        string error = "The Protein does not contain Residue " + c1.str();
        error += insertion;
        throw PANTZ_error (error + "\n");}
    // Return 0 so the function compiles correctly
    return 0;
}

// Get a ResiduePtr object
PROT::ResiduePtr PROT::Protein::operator[] (const size_t i) {
    // Throw an error if the Protein is empty
    if (m_count == 0) {
        string error = "It is not possible to access a Residue of an empty "
                       "Protein.\n";
        throw PANTZ_error (error);}
    // Also throw an error if the index is out of bounds
    else if (i >= m_count) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_count;
        string error = c1.str() + " is not a valid Residue index in a Protein "
                       "of " + c2.str() + " Residues.\n";
        throw PANTZ_error (error);}
    // return the ResiduePtr
    ResiduePtr ptr (&(m_residues[i]));
    return ptr;
}
