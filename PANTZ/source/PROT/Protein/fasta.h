/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the FASTA sequence method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

string PROT::Protein::fasta (const size_t L=0) const {
    // if the protein is empty, throw an error
    if (m_count == 0) {
        string error = "The fasta method does not work for an empty Protein.\n";
        throw PANTZ_error (error);}
    // Store the output here
    string sequence = "";
    // Use a try statement in case any of the residues are not amino acids
    try {
        // A count of how many terms to include in a line
        size_t N = 0;
        // Go through the residues
        for(size_t i=0; i<m_count; ++i) {
            // If needed, add an end line character
            if ((L > 0) && (N >= L)) {sequence += "\n"; N = 0;}
            // Add the 1 letter amino acid code of the residue
            sequence += m_residues[i].AA1();
            ++N;}}
    // If any of the residues were not amino acids
    catch (PANTZ_error& e) {
        string error = "This error occurred in the fasta method of the Protein "
                       "class.\n";
        throw PANTZ_error (e, error);}
    // Now that the sequence has been made, add the header information
    string output = ">Protein ";
    output += m_name;
    output += "\n" + sequence + "\n";
    return output;
}
