/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Structure.h
 * header file, and includes pre-processor directives for that behavior. The
 * Structure class methods for getting content are defined here. */

// Make sure the Structure class is currently being loaded
#ifndef Structure_Loading_Status
#error Structure methods must be included by the Structure.h header file
#endif

// Get a protein
PROT::Protein * PROT::Structure::protein (const size_t i) {
    if (i >= m_proteins.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_proteins.size();
        string error = c1.str() + " is not a valid Protein index in a "
                       "Structures with " + c2.str() + " Proteins.\n";
        throw PANTZ_error (error);}
    return m_proteins[i];
}

// Get a name
string PROT::Structure::name (const size_t i) const {
    if (i >= m_names.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_names.size();
        string error = c1.str() + " is not a valid name index in a "
                       "Structure with " + c2.str() + " names.\n";
        throw PANTZ_error (error);}
    return m_names[i];
}
