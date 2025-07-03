/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the declaration of the ADHOC class. This class is the most
 * complex Protocol in PANTZ. It can implement any Method of the software.
 * Because it follows user-defined instructions, there is not an easy way to
 * error check that all of the instructions make sense. As a result, they're
 * checked as they're run. */

// This file is intended to be included directly from the Protocols.h header
// file
#ifndef Protocol_Loading_Status
#error ADHOC.h must be included in a compiled program by Protocols.h
#endif

// Define the ADHOC class
class PROTOCOL::ADHOC {

    // The information stored in the class
    private:
        // The name of the calculations
        string m_name;
        // Where any output files for these calculations should be stored
        string m_output_path;
        // A summary file of the things the calculations have done
        ofstream m_output;
        // A vector of the PDB files that protein structures are found in
        vector<PROT::PDB> m_pdb_files;
        // A vector of the Protein objects the calculations are working with
        vector<PROT::Protein> m_proteins;
        // An index of the last command in the interface that was actually
        // processed
        size_t m_index;

    // There is a private method used in initializing the class that does a
    // frequent error checking calculation
    private:
        vector<string> * standard_error (const size_t, string&, Interface&);

    // The public class methods
    public:
        // The class destructor makes sure the output file is closed
        ~ADHOC () {if (m_output.is_open()) {m_output.close();}}
        // The class is constructed from an Interface object
        ADHOC (Interface&);
        // Run the calculations
        void run (Interface&);

    // End the class definition
};

// Because the class has few methods, they're implemented directly in this file.

// A standard error checking that occurs repeatedly in the constructor
vector<string> * PROTOCOL::ADHOC::standard_error (const size_t N, 
                                                  string& commandType,
                                                  Interface& interface) {
    // Increment the index
    m_index += 1;
    // Validate that this index is still in the instructions
    if (m_index >= N) {
        string error = "ADHOC PANTZ calculations initialized with no commands.\n";
        throw PANTZ_error (error);}
    // Update the command type and command
    commandType = interface.command_type(m_index);
    return interface.command(m_index);
}

// Implement the class constructor
PROTOCOL::ADHOC::ADHOC (Interface& interface) {
    // Validate that the interface has commands
    size_t N = interface.commands();
    if (N == 0) {
        string error = "Algorithm error: ADHOC calculations should not be "
                       "initialized from an empty interface.\n";
        throw PANTZ_error (error);}
    // Indicate that this is an ADHOC calculation
    string protocolType = "ADHOC";
    // Get the first command
    m_index = 0;
    string commandType = interface.command_type(m_index);
    vector<string> * command = interface.command(m_index);
    // If the command type is wrong, throw an error
    if (commandType != "calculation type") {
        string error = "Algorithm error: the first command for ADHOC "
                       "calculations must be 'calculation type', not '"
                     + commandType + "'\n";
        throw PANTZ_error(error);}
    // Validate that the command has 1 entry
    PROTOCOL::validate_command(command, 1, protocolType, commandType);
    // Convert it to a lower case word and confirm it is adhoc
    string text = (*command)[0]; Text::lower(text);
    if (text != "adhoc") {
        string error = "Algorithm error: ADHOC object being constructed for "
                     + (*command)[0] + " calculations.\n";
        throw PANTZ_error(error);}
    // The second command must be the name of the calculations
    m_index = 1;
    commandType = interface.command_type(m_index);
    command = interface.command(m_index);
    if (commandType != "calculation name") {
        string error = "The second command in ADHOC calculations must be the "
                       "calculation name, not " + commandType + "\n";
        throw PANTZ_error (error);}
    // Validate the command, then store the name
    PROTOCOL::validate_command(command, 1, protocolType, commandType);
    m_name = (*command)[0];
    // Update the index and command inforamtion
    command = standard_error (N, commandType, interface);
    // The next several commands have built in, default behaviors if they're not
    // provided
    if (commandType == "output path") {
        // Validate the command
        PROTOCOL::validate_command(command, 1, protocolType, commandType);
        // Store the information
        m_output_path = (*command)[0];
        // If the path doesn't end with a /, add it
        if (!Text::endswith(m_output_path, '/')) {m_output_path += "/";}
        // Update the index and command information
        command = standard_error (N, commandType, interface);}
    // Otherwise, store the current folder
    else {m_output_path = "./";}
    // Do similar things for whether or not an output folder should be created
    bool make_output;
    if (commandType == "create output folder") {
        PROTOCOL::validate_command (command, 1, protocolType, commandType);
        string value = (*command)[0];
        make_output = PROTOCOL::true_false(value);
        command = standard_error (N, commandType, interface);}
    else {make_output = true;}
    // Whether or not a summary file should be created
    bool make_file;
    if (commandType == "create summary file") {
        PROTOCOL::validate_command (command, 1, protocolType, commandType);
        string value = (*command)[0];
        make_file = PROTOCOL::true_false(value);
        command = standard_error (N, commandType, interface);}
    else {make_file = true;}
    // Whether or not to overwrite previous content when running the
    // calculations
    bool overwrite;
    if (commandType == "overwrite previous") {
        PROTOCOL::validate_command (command, 1, protocolType, commandType);
        string value = (*command)[0];
        overwrite = PROTOCOL::true_false (value);
        command = standard_error (N, commandType, interface);}
    else {overwrite = false;}

    // Validate that the class recognizes all remaining commands
    for(size_t i=m_index; i<N; ++i) {
        commandType = interface.command_type(i);
        if (commandType == "load protein") {continue;}
        else if (commandType == "write proteins") {continue;}
        else if (commandType == "charmm") {continue;}
        else if (commandType == "rosetta") {continue;}
        // if the loop reached this point, that is not a recognized command type
        // for the ADHOC calculations
        string message = "The ADHOC calculations do not recognize the "
                       + commandType + " method.\n";
        throw PANTZ_error (message);}

    // Now actually make folders and files as needed
    if (make_output) {
        // Make the folder
        m_output_path += m_name + "/";
        PROTOCOL::make_folder (m_output_path, overwrite);
        // If output files should also be created, create a copy of the
        // instructions being used for this run.
        if (make_file) {
            // Create a file in it that includes all of the lines of the instruction
            // file
            string fileName = m_output_path + m_name + "_Instructions.txt";
            // Use the output stream to make that file
            PROTOCOL::open_file(fileName, m_output, overwrite);
            for(size_t i=0; i<interface.lines(); ++i) {
                m_output << interface.line(i) << "\n";}
            // Close that output file
            m_output.close();}}
        // Go through the lines of the interface and write them.
    if (make_file) {
        // The name of the output file
        string fileName = m_output_path + m_name + "_Summary.txt";
        // Open it for writing
        PROTOCOL::open_file(fileName, m_output, overwrite);
        // Create a string summarizing the content
        string message = "The " + m_name + " calculations started on "
                       + METHODS::time_stamp () + "\n";
        m_output << message << endl;}
    // End the constructor
}

// The run method of the Interface class carries out the required calculations
void PROTOCOL::ADHOC::run (Interface& interface) {
    // Repeat this for every command
    while (m_index < interface.commands()) {
        // Get the current command
        string commandType = interface.command_type (m_index);
        vector<string> * command = interface.command (m_index);
        // Carry out the command
        if (commandType == "load protein") {
            // Load the proteins
            METHODS::load_protein(command, m_proteins, m_pdb_files, m_output);}
        else if (commandType == "write proteins") {
            // Create the appropriate file location information
            string label = m_output_path + m_name;
            // Call the method
            METHODS::write_proteins(command, m_proteins, label, m_output);}
        else if (commandType == "charmm") {
            // Get how charmm is supposed to be used
            string how = (*command)[0];
            // Convert it to lower case letters
            string use = how; Text::lower(use);
            // run the appropriate CHARMM calculation
            if (use == "fixed backbone energy minimization") {
                CHARMM::Energy_Minimization(m_proteins, m_output, how, 
                                            m_output_path);}
            else if (use == "harmonic backbone energy minimization") {
                CHARMM::Energy_Minimization(m_proteins, m_output, how,
                                            m_output_path);}
            else if (use == "all atom energy minimization") {
                CHARMM::Energy_Minimization(m_proteins, m_output, how,
                                            m_output_path);}
            else if (use == "add missing atoms") {
                CHARMM::Missing_Atoms(m_proteins, m_output, m_output_path);}
            else {
                string error = "The ADHOC calculations do not implement the "
                             + how + " CHARMM method.\n";
                throw PANTZ_error (error);}}
        else if (commandType == "rosetta") {
            string how = (*command)[0];
            string use = how; Text::lower(use);
            if (use == "energy minimization") {
                Rosetta::Energy_Minimization(m_proteins, m_output, 
                                             m_output_path);}
            else if (Text::startswith(use, "interface analysis")) {
                Rosetta::Interface_Analyzer(m_proteins, m_output,
                                            m_output_path, how);}
            else if (use == "per residue") {
                Rosetta::Per_Residue (m_proteins, m_output, m_output_path);}
            else {
                string error = "The ADHOC calculations do not implement the "
                             + how + " Rosetta Method.\n";
                throw PANTZ_error (error);}}
        // If the commadn type wasn't implemented
        else {
            string error = "The ADHOC calculations do not implement the "
                         + commandType + " method.\n";
            throw PANTZ_error (error);}
        // increment the index
        m_index += 1;}
    // When the calculations are done, close the output file if it is open
    if (m_output.is_open()) {
        string message = "Calculations ended on " + METHODS::time_stamp() 
                       + "\n";
        m_output << message;
        m_output.close();}
}
