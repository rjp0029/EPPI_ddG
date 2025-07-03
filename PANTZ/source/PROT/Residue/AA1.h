/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * AA1 method of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// The AA1 function of the Residue class
char PROT::Residue::AA1 () const {
    // Confirm that the Residue is an amino acid
    if (!is_amino_acid()) {
        string error = "The Residue AA1 method only works for amino acids, not: "
                     + m_name + "\n";
        throw PANTZ_error (error);}
    // Get it's name
    string label = m_name;
    // If the name is an alternative histidine label, change it
    if ((label == "HSD") || (label == "HSE")) {label = "HIS";}
    // Find the index of the amino acid
    for(size_t i=0; i<20; ++i) {
        // If the name matches
        if (PROT::AA3[i] == label) {return PROT::AA1[i][0];}}
    // If the function reaches this point, something went very wrong
    string error = "Major algorithm error regarding amino acid names.\n";
    throw PANTZ_error (error);
    // Return 'Z' as a placeholder
    return 'Z';
}
