/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the construction method of the EPPI::Feature class */

// Make sure it is being loaded from the expected file
#ifndef EPPI_Feature_Loading_Status
#error EPPI::Feature::constructor must be included by Feature.h
#endif

// Implement the function
EPPI::Feature::Feature (const string& label, const size_t N) {
    // Set the name of the feature, then strip whitespace and make sure it is
    // lowercase
    m_name = label; Text::strip(m_name); Text::lower(m_name);
    // And the index
    m_index = N;
    // Separate the name into primary feature labels
    string part = "";
    // Go through the name's characters
    for(size_t i=0; i<m_name.size(); ++i) {
        char L = m_name[i];
        // If it is a math operation
        if (((L == '*') || (L == '/')) || ((L == '+') || (L == '-'))) {
            // If no primary feature has been identified, throw an error
            if (part.size() == 0) {
                string error = "Invalid Feature name: " + label + "\n";
                throw PANTZ_error (error);}
            // Otherwise
            else {
                // Store the part as a primary feature name - the validity of
                // that will be checked later
                m_primary.push_back(part);
                // Reset the part to an empty string
                part = "";
                // Store the operation
                m_operations.push_back(L);}}
        // If it isn't a math operation, add it to the primary feature
        else {part += L;}}
    // Store the last part in the primary features list
    if (part.size() > 0) {m_primary.push_back(part);}
    else {
        string error = "Invalid Feature name: " + label + "\n";
        throw PANTZ_error (error);}
    // If there are no primary parts, throw an error
    if (m_primary.size() == 0) {
        string error = "Feature construction with an invalid name\n";
        throw PANTZ_error (error);}
    // Go through each primary part
    for(size_t i=0; i<m_primary.size(); ++i) {
        // If it is a known primary feature
        bool known = false;
        // Loop through the primary features
        for(size_t j=0; j<EPPI::PrimaryFeaturesList.size(); ++j) {
            // If the name matches
            if (m_primary[i] == EPPI::PrimaryFeaturesList[j]) {
                known = true;
                break;}}
        // If it isn't known, throw an error
        if (!known) {
            string error = "Unknown primary feature (" + m_primary[i]
                         + ") encountered during Feature construction\n";
            throw PANTZ_error (error);}}
}
