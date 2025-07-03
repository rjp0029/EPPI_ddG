/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the is_hydrogen method of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// Determine whether or not an Atom is a hydrogen
bool PROT::Atom::is_hydrogen () const {
    // Loop through the characters in the Atom's name
    if (m_name.size() > 0) {
        for(size_t i=0; i<m_name.size(); ++i) {
            // If it is a digit, continue
            if (Text::is_digit(m_name[i])) {continue;}
            // For the function to have reached this point it is looking at the
            // first non-digit character in the atom's name.
            return (m_name[i] == 'H');}}
    // If the function reached this point, there was no letter char in the name,
    // so it is not a hydrogen
    return false;
}

// Determine an Atom's chemical element
char PROT::Atom::determine_element () const {
    // If the name is empty, throw an error
    if (m_name.size() == 0) {
        string error = "The determine element method of the Atom class does "
                       "not work for an unnamed Atom.\n";
        throw PANTZ_error (error);}
    // Go through the characters of the atom's name
    for (size_t i=0; i<m_name.size(); ++i) {
        // If it is a digit, continue
        if (Text::is_digit(m_name[i])) {continue;}
        // Otherwise, return the character
        return m_name[i];}
    // If the function reached this point, throw an error
    string error = "The determine element method of the Atom class failed for "
                 + m_name + "\n";
    throw PANTZ_error (error);
    // Return a variable so the function compiles
    return 'X';
}
