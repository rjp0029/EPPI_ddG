/* Created by the Pantazes Lab at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the standard construction method of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// The standard constructor of the Atom class uses a string of text
PROT::Atom::Atom (const string& input) {
    // Assign initial values to the Atom
    initialize ();
    // Create a copy of the string that can be modified and analyzed
    string line = input;
    // Strip the whitespace and convert it to all capital letters
    Text::strip(line); Text::upper(line);
    // Put everything in a try statement
    try {
        // There must be at least 66 columns
        if (line.size () < 66) {
            string error = "This text is too short to contain a PDB atom.\n";
            throw PANTZ_error (error);}
        // Determine the Atom's type
        if (Text::startswith(line, "ATOM  ")) {m_type = "ATOM";}
        else if (Text::startswith(line, "HETATM")) {m_type = "HETATM";}
        else {
            string error = "PDB atom lines must start with ATOM or HETATM.\n";
            throw PANTZ_error (error);}
        // Get the Atom's number
        string value = line.substr(6, 5); Text::strip(value);
        if (Text::is_integer(value)) {
            stringstream c; c << value; c >> m_number;
            CHECK::atom_number (m_number);}
        else {
            string error = "Columns 7-11 do not contain an integer.\n";
            throw PANTZ_error (error);}
        // Get the Atom's name
        m_name = line.substr(12, 4); Text::strip(m_name);
        CHECK::atom_name(m_name);
        // Get the alternative location information
        m_alt = line[16];
        CHECK::alt_location(m_alt);
        // Get the residue's name
        m_residue = line.substr(17, 3); Text::strip(m_residue);
        CHECK::residue_name (m_residue);
        // Get the Protein's name
        m_protein = line[21];
        CHECK::protein_name(m_protein);
        // Get the Residue's number
        value = line.substr(22, 4); Text::strip(value);
        if (Text::is_integer(value)) {
            stringstream c; c << value; c >> m_residue_number;
            CHECK::residue_number (m_residue_number);}
        else {
            string error = "Columns 23-26 do not contain an integer.\n";
            throw PANTZ_error (error);}
        // Get the Residue's insertion code
        m_insertion = line[26];
        CHECK::insertion_code (m_insertion);
        // Get the Atom's coordinates (note that AtomCoordinates is not used
        // here because PDB files ONLY have 3 coordinates)
        for(size_t i=0; i<3; ++i) {
            value = line.substr(30 + (i*8), 8); Text::strip(value);
            if (Text::is_number(value)) {
                stringstream c; c << value; c >> m_coors[i];
                CHECK::atom_coordinate(m_coors[i]);}
            else {
                string error = "PDB atom coordinates have to be numbers.\n";
                throw PANTZ_error (error);}}
        // Get and store the Atom's occupancy
        value = line.substr(54, 6); Text::strip(value);
        if (Text::is_number(value)) {
            stringstream c; c << value; c >> m_occupancy;
            CHECK::occupancy(m_occupancy);}
        else {
            string error = "Columns 55-60 do not contain a number.\n";
            throw PANTZ_error (error);}
        // Get and store the Atom's temperature
        value = line.substr(60, 6); Text::strip(value);
        if (Text::is_number(value)) {
            stringstream c; c << value; c >> m_temperature;
            CHECK::temperature(m_temperature);}
        else {
            string error = "Columns 61-66 do not contain a number.\n";
            throw PANTZ_error (error);}
        // The element and charge information may not be present
        if (line.size() >= 78) {
            m_element = line.substr(76, 2); Text::strip(m_element);
            CHECK::element(m_element);
            if (line.size() >= 80) {
                m_charge = line.substr(78, 2); Text::strip(m_charge);
                CHECK::charge (m_charge);}}
    // Catch and handle any errors that occurred during this process
    } catch (PANTZ_error& e) {
        // Store the initial input and strip it of whitespace
        string data = input; Text::strip(data);
        // Throw a new error with the same message but also the input line
        throw PANTZ_error (e, data);}
    // End the function
}
