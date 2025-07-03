/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the check file status method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Check to see whether or not a PDB file is obsolete or theoretical
void PROT::PDB::check_file_status () {
    // Set the default values to false
    m_obsolete = false;
    m_theoretical = false;
    // Search through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line indicates that the file is obsolete
        if (Text::startswith(m_lines[i], "OBSLTE")) {
            m_obsolete = true;}
        // If the line is a REMARK or ATOM line, the search can be finished
        // because OBSLTE entries are required to come BEFORE that information
        else if ((Text::startswith(m_lines[i], "REMARK")) ||
                 (Text::startswith(m_lines[i], "ATOM"))) {break;}}
    // Also check to see if the file is theoretical
    m_theoretical = Text::contains(m_type, "THEORETICAL");
}
