/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the copy method of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// Copy the information from another instance of the PDB class
void PROT::PDB::copy (const PDB * other) {
    // Copy the non-vector variables
    m_name = other->m_name;
    m_folder = other->m_folder;
    m_type = other->m_type;
    m_resolution = other->m_resolution;
    m_obsolete = other->m_obsolete;
    m_theoretical = other->m_theoretical;
    // Copy the lines vector
    m_lines.clear(); m_lines.reserve(other->m_lines.size());
    if (other->m_lines.size() > 0) {
        for(size_t i=0; i<other->m_lines.size(); ++i) {
            m_lines.push_back(other->m_lines[i]);}}
    // Copy the proteins
    m_proteins.clear(); m_proteins.reserve(other->m_proteins.size());
    if (other->m_proteins.size() > 0) {
        for(size_t i=0; i<other->m_proteins.size(); ++i) {
            m_proteins.push_back(other->m_proteins[i]);}}
    // Copy the Structures, making sure the pointers will point to the
    // Proteins in this PDB file
    m_structures.clear(); m_structures.reserve(other->m_structures.size());
    if (other->m_structures.size() > 0) {
        for(size_t i=0; i<other->m_structures.size(); ++i) {
            m_structures.push_back(Structure(other->m_structures[i], m_proteins));}}
}
