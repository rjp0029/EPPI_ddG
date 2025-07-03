/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the atom collection method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Collect the Atoms that make up the Proteins
void PROT::PDB::collect_Atoms (vector<vector<Atom> >& atoms) {
    // Set up the atoms vector to be accessible by character indexing
    atoms.clear(); atoms.resize(128);
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with ATOM or HETATM
        if ((Text::startswith(m_lines[i], "ATOM")) ||
            (Text::startswith(m_lines[i], "HETATM"))) {
            // Try to create an atom from the line
            try {
                Atom atom (m_lines[i]);
                // Convert the Atom's character to an integer
                int n = (int) atom.protein();
                // If the Atom's alternative location characteristic is 'A' or
                // ' ', store it
                if (atom.m_alt == 'A') {
                    atom.m_alt = ' '; atoms[n].push_back(atom);}
                else if (atom.m_alt == ' ') {
                    atoms[n].push_back(atom);}}
            // If an error occurs, just ignore it
            catch (PANTZ_error) {continue;}}}
}
