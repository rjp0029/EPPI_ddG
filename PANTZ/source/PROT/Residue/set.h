/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * methods for setting number / chain of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Set the numbering information of the Residue without error checking
void PROT::Residue::private_set_number (const long n, const char i = ' ', 
                                        const bool internal = true) {
    // If this is a modification of the Residue's internal number information,
    // then only the internal number is set
    if (internal) {m_internal = n;}
    // Otherwise, update the number and insertion information
    else {m_number = n; m_insertion = i;}
    // Note that this function does NOT change the information in Atoms. That is
    // done when the Atoms are being output for something.
}

// Set the number with error checking
void PROT::Residue::set_number (const long n, const char i = ' ',
                                const bool internal = true) {
    // Check the provided information
    CHECK::residue_number (n);
    CHECK::insertion_code (i);
    // Set the information
    private_set_number (n, i, internal);
}

// Set the protein information of the Residue without error checking
void PROT::Residue::private_set_protein (const char L) {
    // Store the value of the Residue
    m_protein = L;
    // If there are Atoms, update them, too
    if(m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_atoms[i].m_protein = L;}}
}

// Set the protein information wiht error checking
void PROT::Residue::set_protein (const char L) {
    CHECK::protein_name (L);
    private_set_protein(L);
}

// The set number function of the ResiduePtr class is declared here to make sure
// that the default behaviors match
void PROT::ResiduePtr::set_number (const long n, const char i = ' ',
                                   const bool internal = true) {
    // Check the residue ptr
    check ();
    // Set it's number
    m_ptr->set_number(n, i, internal);
}
