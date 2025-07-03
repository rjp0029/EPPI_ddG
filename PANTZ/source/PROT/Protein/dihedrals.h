/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the calculate_dihedrals method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Calculate the dihedral angles of the Protein's residues
void PROT::Protein::calculate_dihedrals() {
    // If the protein is empty, throw an error
    if (m_count == 0) {
        string error = "The calculate_dihedrals method does not work for an "
                       "empty Protein.\n";
        throw PANTZ_error (error);}
    // Confirm every residue is an amino acid
    for (size_t i=0; i<m_count; ++i) {
        if (!m_residues[i].is_amino_acid()) {
            string error = "The calculate_dihedrals method only works when "
                           "all Residues in a Protein are amino acids.\n";
            throw PANTZ_error (error);}}
    // Specify strings of the backbone atoms
    string name1 = "N";
    string name2 = "CA";
    string name3 = "C";
    // Use a try statement in case any residue is missing an atom
    try {
        // Loop through the Protein's Residues
        for(size_t i=0; i<m_count; ++i) {
            // Get the 3 backbone Atoms
            AtomPtr N = m_residues[i][name1];
            AtomPtr CA = m_residues[i][name2];
            AtomPtr C = m_residues[i][name3];
            // If this is not the first Residue, calculate the phi dihedral 
            // angle using the C from the previous residue
            if (i > 0) {
                AtomPtr Cp = m_residues[i-1][name3];
                m_residues[i].m_phi = calculate_dihedral(Cp, N, CA, C);}
            // If this is not the last residue
            if (i < m_count - 1) {
                // Get the Nitrogen in the next Residue
                AtomPtr Nn = m_residues[i+1][name1];
                // Calculate psi
                m_residues[i].m_psi = calculate_dihedral(N, CA, C, Nn);
                // Get the alpha carbon of the next residue
                AtomPtr CAn = m_residues[i+1][name2];
                // Calculate omega
                m_residues[i].m_omega = calculate_dihedral(CA, C, Nn, CAn);}}}
    // Catch errors 
    catch (PANTZ_error& e) {
        string error = "This error occurred in the calculate_dihedrals method "
                       "of the Protein class.\n";
        throw PANTZ_error (e, error);}
}

