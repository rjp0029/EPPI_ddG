/* Created by the Pantazes Lab at Auburn University.
 *
 * This file runs CHARMM energy minimizations in the PANTZ software suite. */

// It must be included by CHARMM.h
#ifndef CHARMM_Loading_Status
#error CHARMM methods must be included by CHARMM.h
#endif

// Implement the method to run energy minimizations
void CHARMM::Energy_Minimization (vector<PROT::Protein>& proteins,
             ofstream& output, const string& how, const string& path) {
    // Store the script here
    string script = "* ";
    // A message that will be part of the script and possibly output to the
    // summary file
    string message = "";
    // Whether harmonic or fixed backbone calculations are happening
    bool harmonic = false;
    bool fixed = false;
    // Go through the options of how to do the calculations
    string use = how;
    Text::lower(use);
    if (use == "fixed backbone energy minimization") {
        fixed = true;
        message = "CHARMM Fixed Backbone Energy Minimization";}
    else if (use == "harmonic backbone energy minimization") {
        harmonic = true;
        message = "CHARMM Harmonic Backbone Energy Minimization";}
    else if (use == "all atom energy minimization") {
        message = "CHARMM All Atom Energy Minimization";}
    else {
        string error = "The CHARMM::Energy_Minimization function does not "
                       "recognize the " + how + " method.\n";
        throw PANTZ_error(error);}
    // Update the script
    script += message + "\n";
    // If the summary file is open, write to it
    if (output.is_open()) {
        output << message << " started on " << METHODS::time_stamp() << endl;}
    // Get the current working directory
    string cwd = get_cwd(path);
    // Change to the output path folder
    change_directory(path);
    // Make the rest of the script
    warning_bomb (script);
    load_topology_parameter(script);
    proteins_for_charmm(script, proteins);
    add_missing_atoms (script);
    energy_minimization(script, harmonic, fixed);
    output_proteins (script, proteins);
    script += "stop\n";
    // Run the script
    string label = "energy_minimization";
    run_charmm_script(script, label);
    // Load the proteins
    proteins_from_charmm(proteins);
    // Clean up the files
    clean_up (label, proteins, true, true);
    // Change back to the original directory
    change_directory(cwd);
    // if the output file is open, update it
    if (output.is_open()) {
        output << "Calculations ended on " << METHODS::time_stamp() << "\n" << endl;}
}
