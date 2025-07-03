/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the exclusions method of the Atom class. This method is deprecated, but is
 * being kept in the code base in case it is ever needed in the future. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Determine whether or not non-bonded exclusions should be applied in the
// calculation of energies between 2 Atoms
int PROT::Atom::exclusions (const Atom * other) const {
    // As currently implemented, this function is ONLY correct if THIS atom is
    // a side-chain atom in a Rotamer. Currently, this is the only instance
    // when intra-protein energies are calculated. If that is ever changed,
    // this function will need to be changed significantly

    // If the two Atoms are in different proteins, they definitely don't need
    // non-bonded exclusions
    if (m_protein != other->m_protein) {return 0;}
    // If the two Atoms are at least 2 Residues apart from one another,
    // non-bonded exclusions aren't needed
    else if ((other->m_residue_number < m_residue_number - 1) || 
             (other->m_residue_number > m_residue_number + 1)) {return 0;}
    // If the two atoms are in the same residue and are both side chain atoms,
    // make sure energy calculations are NOT done
    else if ((other->m_residue_number == m_residue_number) && 
            ((!other->is_backbone_atom()) && (!is_backbone_atom()))) {return -1;}
    // If the other atom is in the previous residue
    else if (other->m_residue_number == m_residue_number - 1) {
        // Proline gets special handling, because its delta carbon is bound to
        // its nitrogen. This opens up the possibility of non-bonded
        // exclusions that aren't available to other amino acids
        if (m_residue == "PRO") {
            // If this is the delta carbon
            if (m_name == "CD") {
                // If the other atom is C
                if (other->m_name == "C") {return 3;}
                else if ((other->m_name == "CA") || (other->m_name == "O")) {
                    return 4;}}
            // If this is the Gamma carbon or a delta carbon hydrogen
            else if ((m_name == "CG") || ((m_name == "HD1") || 
                     (m_name == "HD2"))) {
                // In that case, a non-bonded exclusion of 4 is used with the
                // C from the previous residue
                if (other->m_name == "C") {return 4;}}}
        // Even for prolines, if this atom is the beta carbon (or HA1 for
        // glycines) and the other Atom is C, a non-bonded exclusion of 4 is
        // correct
        if (((m_name == "CB") || (m_name == "HA1")) && (other->m_name == "C")) {
            return 4;}
        // In all other cases where the other atom is in the previous residue,
        // non-bonded exclusions are not appropriate
        return 0;}
    // If the other atom is in the next residue
    else if (other->m_residue_number == m_residue_number + 1) {
        // If this atom is the beta carbon or HA1 and the next atom is a
        // nitrogen, a non-bonded exclusion of 4 is correct
        if (((m_name == "CB") || (m_name == "HA1")) && (other->m_name == "N")) {
            return 4;}
        // In all other cases, no non-bonded exclusion is appropriate
        return 0;}
    // If the function has reached this point, we know that the the two Atoms
    // are in the same residue in the same protein. Because of the delta
    // carbon to nitrogen bond of prolines, they need special handling
    else if (m_residue == "PRO") {
        // Non-proline, N-terminal amino acids that are having proline
        // considered as a rotamer have an "extra" atom. Use the smallest
        // non-bonded exclusion value without question to ensure no energies
        // are calculated with that atom
        if (other->m_name == "HT3") {return 2;}
        // A subsequent portion of this function will calculate the values
        // from atoms off of the alpha carbon. Here, do the calculations for
        // atoms off of the nitrogen
        else if (other->m_name == "N") {
            if (m_name == "CD") {return 2;}
            else if ((m_name == "CG") || 
                    ((m_name == "HD1") || (m_name == "HD2"))) {return 3;}
            // The beta carbon is a 3, going through the alpha carbon. That
            // gets picked up later.
            else if ((m_name == "HG1") || (m_name == "HG2")) {return 4;}}
        // Backbone Atoms bound to the nitrogen
        else if ((((other->m_name == "CA") || (other->m_name == "HN")) ||
                  ((other->m_name == "HN1") || (other->m_name == "HN2"))) ||
                  ((other->m_name == "HT1") || (other->m_name == "HT2"))) {
            if (m_name == "CD") {return 3;}
            // The gamma carbon with the alpha carbon is a 3 via the beta
            // carbon. Get that correct here
            else if ((m_name == "CG") || (other->m_name == "CA")) {return 3;}
            else if ((m_name == "CG") || 
                    ((m_name == "HD1") || (m_name == "HD2"))) {return 4;}}
        // Backbone Atoms bound to those atoms in this residue
        else if ((other->m_name == "C") || 
                ((other->m_name == "HA") || (other->m_name == "HA2"))) {
            if (m_name == "CD") {return 4;}}}
    // Search through the possibilities for within the same residue via
    // alpha-carbon connections. This is also done for prolines, since they
    // still have the same alpha-carbon connections
    // Rather than doing many complicated if statements, arrays of values will
    // be checked for matching residue names. First are the backbone atom
    // names
    string first [1] = {"CA"};
    string second [4] = {"C", "N", "HA", "HA2"};
    string third [9] = {"O", "OT1", "OT2", "HN", "HN1", "HN2", "HT1", "HT2",
                        "HT3"};
    // And now the side chain atom names
    string beta [2] = {"CB", "HA1"};
    string gamma [10] = {"CG", "CG1", "CG2", "HB", "HB1", "HB2", "HB3", "OG", 
                         "OG1", "SG"};
    string delta [17] = {"CD", "CD1", "CD2", "HG", "HG1", "HG2", "HG11", "HG12",
                         "HG13", "HG21", "HG22", "HG23", "OD1", "OD2", "ND1",
                         "ND2", "SD"};
    // Check the Beta Carbon Atoms
    for(size_t i=0; i<2; ++i) {
        if (m_name == beta[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 2;}}
            for(size_t j=0; j<4; ++j) {
                if (other->m_name == second[j]) {return 3;}}
            for(size_t j=0; j<9; ++j) {
                if (other->m_name == third[j]) {return 4;}}}}
    // Gamma position atoms
    for(size_t i=0; i<10; ++i) {
        if (m_name == gamma[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 3;}}
            for(size_t j=0; j<4; ++j) {
                if (other->m_name == second[j]) {return 4;}}}}
    // Delta position atoms
    for(size_t i=0; i<17; ++i) {
        if (m_name == delta[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 4;}}}}
    // In any other circumstance, return 0
    return 0;
}
