/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the str method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Create a formatted string of text containing the Protein's contents
string PROT::Protein::str (const bool internal = false) {
    // STore the output here
    string output; output.reserve(number_of_atoms() * AtomStringLength);
    // If there are Residues
    if (m_count > 0) {
        // Loop through the Residues
        for(size_t i=0; i<m_count; ++i) {
            // Add the Residue's information to the Protein string
            output.append(m_residues[i].str(internal));}}
    return output;
}
