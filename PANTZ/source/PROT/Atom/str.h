/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the str method of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Generate a string of PDB-formatted text containing the Atom's information
string PROT::Atom::str () const {
    // The output is stored here
    string output;
    // Reserve the appropriate number of characters
    output.reserve(AtomStringLength);
    // Include the Atom's attributes with proper formatting and spacing
    // The Atom's type
    Text::ljust_insert (output, m_type, 6, ' ');
    // The Atom's number
    stringstream c1; c1 << m_number;
    Text::rjust_insert(output, c1.str(), 5, ' ');
    // A blank space
    output.push_back(' ');
    // The Atom's name
    if (m_name.size() < 4) {
        output.push_back(' '); Text::ljust_insert(output, m_name, 3, ' ');}
    else {Text::ljust_insert(output, m_name, 4, ' ');}
    // The alternate location character
    output.push_back(m_alt);
    // The Residue's name
    Text::ljust_insert(output, m_residue, 3, ' ');
    // A blank space
    output.push_back (' ');
    // The protein's name
    output.push_back(m_protein);
    // The Residue's number
    stringstream c2; c2 << m_residue_number;
    Text::rjust_insert(output, c2.str(), 4, ' ');
    // The residue's insertion code
    output.push_back(m_insertion);
    // Three blank spaces
    for(size_t i=0; i<3; ++i) {output.push_back(' ');}
    // The Atom's coordinates. As with loading from the PDB file, 3
    // coordinates are explicitly used here
    for(size_t i=0; i<3; ++i) {
        stringstream c3; c3 << fixed << setprecision(3) << m_coors[i];
        Text::rjust_insert(output, c3.str(), 8, ' ');}
    // Occupancy
    stringstream c4; c4 << fixed << setprecision(2) << m_occupancy;
    Text::rjust_insert(output, c4.str(), 6, ' ');
    // Temperature
    stringstream c5; c5 << fixed << setprecision(2) << m_temperature;
    Text::rjust_insert(output, c5.str(), 6, ' ');
    // The Atom's element
    Text::rjust_insert(output, m_element, 12, ' ');
    // The Atom's charge
    Text::rjust_insert(output, m_charge, 2, ' ');
    // Add an end line character to terminate the string
    output.push_back('\n');
    return output;
}

// The rosetta string function reformats an Atom's information to make it work
// in Rosetta
string PROT::Atom::rosetta_str () const {
    // Store the output here
    string output = "";
    // If the atom is a hydrogen, don't do anything else
    if (is_hydrogen()) {return output;}
    // First, set the output equal to the standard str output
    output = str();
    // Reformat it
    output = output.substr(0, 66) + "           ";
    // Add the atom's element
    output += determine_element();
    // End the line
    output += " \n";
    return output;
}
