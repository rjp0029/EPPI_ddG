/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements constructors of the BCProps class that take data
 * formatted from the class's good and row string methods */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error str_constructors.h must be included by BCProps.h
#endif

// The method to load data from the good string function
void EPPI::BCProps::good_str_constructor (const vector<string>& lines) {
    // Initialize the class
    initialize();
    // Validate that there are the expected number of lines
    if (lines.size() != BaseFeaturesList.size() + 1) {
        string error = "BCProps class being initialized from an improperly "
                       "sized vector of strings\n";
        throw PANTZ_error (error);}
    // Allocate memory for the base features
    m_baseFeatures.resize(BaseFeaturesList.size());
    // Make a vector for the base features that are found
    vector<bool> found; found.resize(BaseFeaturesList.size());
    for(size_t i=0; i<BaseFeaturesList.size(); ++i) {found[i] = false;}
    // Loop through the lines
    for(size_t i=0; i<lines.size(); ++i) {
        // Copy the line, strip whitespace, and split it into parts
        string text = lines[i]; Text::strip(text); 
        vector<string> parts; Text::split(parts, text);
        // If this is the first line, it should be the name of the instance
        if (i == 0) {
            // There should be only one part
            if (parts.size() == 1) {m_name = parts[0];}
            else {
                string error = "BCProps instance being constructed from an "
                               "invalid set of strings\n";
                throw PANTZ_error (error);}}
        // Otherwise it should be a base feature and value
        else {
            // Make sure there are two parts
            if (parts.size() != 2) {
                string error = "BCProps instance being constructed from an "
                               "invalid set of strings\n";
                throw PANTZ_error (error);}
            // Get the label of the base feature and make sure it is lower case
            string label = parts[0]; Text::lower(label);
            // Get the index of that base feature
            size_t I = base_feature_index (label);
            // If that base feature has already been set, throw an error
            if (found[I] == true) {
                string error = "BCProps instance being constructed from a "
                               "list where " + label + " is listed multiple "
                               "times\n";
                throw PANTZ_error (error);}
            // Mark the property as found
            found[I] = true;
            // Get the value and store it.
            stringstream c1; c1 << parts[1]; c1 >> m_baseFeatures[I];}}
    // End the class
}

// The method to load data from the row string output
void EPPI::BCProps::row_str_constructor (const string& label, 
           const vector<string>& parts) {
    // Initialize the class
    initialize ();
    // Validate the number of entries
    if (parts.size() != BaseFeaturesList.size() + 1) {
        string error = "BCProps class being initialized from an improperly "
                       "sized vector of strings\n";
        throw PANTZ_error (error);}
    // Allocate memory for the base features
    m_baseFeatures.resize(BaseFeaturesList.size());
    // Set the name of the instance
    m_name = label;
    // Make sure it matches the first part
    if (parts[0] != m_name) {
        string error = "BCProps instance being constructed from an invalid "
                       "vector of strings\n";
        throw PANTZ_error (error);}
    // Store the values
    for(size_t i=0; i<BaseFeaturesList.size(); ++i) {
        // Store the numerical value
        stringstream c1; c1 << parts[i+1]; c1 >> m_baseFeatures[i];}
}
