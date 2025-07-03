/* Created by the Pantazes Lab at Auburn University.
 *
 * This function uses CHARMM to add missing atoms to proteins. */

// It must be included by CHARMM.h
#ifndef CHARMM_Loading_Status
#error CHARMM methods must be included by CHARMM.h
#endif

// Implement the method to run energy minimizations
void CHARMM::Missing_Atoms (vector<PROT::Protein>& proteins,
                            ofstream& output, const string& path) {
    // If the output file is open, add a message to it
    string message;
    if (output.is_open()) {
        message = "Adding missing atoms with CHARMM";
        output << message << " started on " << METHODS::time_stamp() << endl;}
    // Get the current working directory
    string cwd = get_cwd(path);
    // Change to the folder for the calculations
    change_directory(path);
    // Store the script here
    string script = "* Add Missing Atoms\n";
    // Make the rest of the script
    warning_bomb (script);
    load_topology_parameter(script);
    proteins_for_charmm(script, proteins);
    add_missing_atoms (script);
    output_proteins (script, proteins);
    script += "stop\n";
    // Run the script
    string label = "missing_atoms";
    run_charmm_script(script, label);
    // Load the proteins
    proteins_from_charmm(proteins);
    // Clean up the files
    clean_up (label, proteins, true, true);
    // Move back to the original folder
    change_directory(cwd);
    // if the output file is open, update it
    if (output.is_open()) {
        output << "Calculations ended on " << METHODS::time_stamp() << "\n" << endl;}
}
