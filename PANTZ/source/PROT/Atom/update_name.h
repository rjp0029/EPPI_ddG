/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the Rosetta naming update methods of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Update an Atom's name for use in Rosetta
void PROT::Atom::update_name_for_Rosetta (bool lastResidue = false) {
    // Because hydrogens are never output to Rosetta, their names don't have to
    // be updated for it. Only update heavy atoms
    if ((m_residue == "ILE") && (m_name == "CD")) {m_name = "CD1";}
    else if ((lastResidue) && (m_name == "OT1")) {m_name = "O";}
    else if ((lastResidue) && (m_name == "OT2")) {m_name = "OXT";}
}

// Update an Atom's name after Rosetta. This is much more complicated because
// Rosetta uses a different Hydrogen naming convention than CHARMM or the
// Rotamer library
void PROT::Atom::update_name_after_Rosetta (bool lastResidue = false) {
    // Update the delta carbon of isoleucine
    if ((m_residue == "ILE") && (m_name == "CD1")) {m_name = "CD";}
    // Terminal oxygens
    else if ((lastResidue) && (m_name == "O")) {m_name = "OT1";}
    else if ((lastResidue) && (m_name == "OXT")) {m_name = "OT2";}
    // Hydrogens are more complicated. While the general rule is to move the
    // leading digit to the end of the name, there are other name schemes, too
    else if (is_hydrogen()) {
        // N-terminal hydrogens
        if ((m_name == "1H") || ((m_name == "2H") || (m_name == "3H"))) {
            char digit = m_name[0];
            m_name = "HT";
            m_name += digit;}
        // Hydrogens on N
        else if (m_name == "H") {m_name = "HN";}
        // If the name starts with a digit, move it to the end of the name
        else if (Text::is_digit(m_name[0])) {
            // Get the digit
            char digit = m_name[0];
            // The number of characters to include in the moved piece
            size_t n = m_name.size() - 1;
            // Get the last n characters of the atom's name
            m_name = m_name.substr(1, n);
            // Add the digit to the end of the name
            m_name += digit;}}
}
