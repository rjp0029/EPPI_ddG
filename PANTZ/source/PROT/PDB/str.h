/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the str(ing) method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// A string representation of the PDB file's information
string PROT::PDB::str () const {
    // Store the string here
    string output = "PDB File: " + m_name + "\n";
    // If it is obsolete, write that and be done
    if (m_obsolete) {
        output += "THIS FILE IS OBSOLETE AND NOT APPROPRIATE FOR USE.\n";
        return output;}
    // List the experiment type
    output += "Experiment Type: " + m_type + "\n";
    // List the resolution
    output += "Resolution: ";
    if (m_resolution > 0) {
        stringstream c1; c1 << fixed << setprecision(3) << m_resolution;
        output += c1.str() + " Angstroms.\n";}
    else {output += " Not Applicable.\n";}
    // List the structures in the PDB file
    output += "\n";
    for(size_t i=0; i<m_structures.size(); ++i) {
        output += m_structures[i].str();}
    return output;
}
