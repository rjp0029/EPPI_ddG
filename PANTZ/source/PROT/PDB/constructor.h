/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the PDB.h header file. It implements
 * the constructor of the PDB class. */

// Confirm that the PDB class is loading the content
#ifndef PDB_Loading_Status
#error PDB methods must be included by PDB.h
#endif

// The constructor of the PDB class
PROT::PDB::PDB (const string& fileName, const string path = "") {
    // Set the name
    m_name = fileName;
    // Set the folder
    m_folder = path;
    // Load the contents of the file
    load();
    // Use a try statement to add to any error that comes up in the subsequent
    // steps
    try {
        identify_type();
        check_file_status();
        if (m_obsolete) {return;}
        identify_resolution ();
        construct_Proteins ();
        create_Structures ();
        }
    catch (PANTZ_error& e) {
        string error = "This occurred in PDB file " + m_name + "\n";
        throw PANTZ_error (e, error);}
}
