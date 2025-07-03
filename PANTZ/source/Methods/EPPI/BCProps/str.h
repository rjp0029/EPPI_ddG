/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the string methods of the BCProps class */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error str.h must be included by BCProps.h
#endif

// Implement the "good string" method
string EPPI::BCProps::good_str () const {
    // Store the output here
    string output = "Binding Complex Name: " + m_name + "\n";
    // List each base feature
    for(size_t i=0; i<BaseFeaturesList.size(); ++i) {
        // Put the value in a stringstream
        stringstream c1; c1 << fixed << setprecision(4) << m_baseFeatures[i];
        // Make the label
        string label = BaseFeaturesList[i] + " ";
        // Insert it in the string
        Text::ljust_insert(output, label, 15, ' ');
        // Add the number
        output += c1.str() + "\n";}
    // Return the string
    return output;
}

// Implement the "row string" method
string EPPI::BCProps::row_str () const {
    // Store the output here
    string output = m_name;
    // Do this for each feature
    for(size_t i=0; i<m_baseFeatures.size(); ++i) {
        // Put the value in a string stream with proper precision
        stringstream c1; c1 << fixed << setprecision(4) << m_baseFeatures[i];
        // Add it to the string
        output += " " + c1.str();}
    output += "\n";
    return output;
}
