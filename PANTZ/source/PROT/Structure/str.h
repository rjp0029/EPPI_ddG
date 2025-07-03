/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Structure.h
 * header file, and includes pre-processor directives for that behavior. The
 * str method of the Structure class are defined here. */

// Make sure the Structure class is currently being loaded
#ifndef Structure_Loading_Status
#error Structure methods must be included by the Structure.h header file
#endif

// Summarize the Structure
string PROT::Structure::str () const {
    // Store the output here
    string output = "Molecule\nName(s):";
    // if there aren't any
    if (m_names.size() == 0) {output += " Unknown\n";}
    // List the names
    else {
        for(size_t i=0; i<m_names.size(); ++i) {
            output += " " + m_names[i];
            if (i < m_names.size() - 1) {output += ";";}
            else {output += "\n";}}}
    // If there are no proteins
    if (m_proteins.size() == 0) {output += "Chains: None\n";}
    // Otherwise, list them
    else {
        for(size_t i=0; i<m_proteins.size(); ++i) {
            stringstream c1; 
            c1 << fixed << setprecision(4) << m_proteins[i]->score();
            output += "Chain ";
            output += m_proteins[i]->name();
            output += ": Score = " + c1.str() + "\n";}}
    output += "\n";
    return output;
}
