/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements methods of the BCProps class that use base and primary
 * features to calculate the primary and derived features */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error calculate_features.h must be included by BCProps.h
#endif

// Implement the primary features calculation method
void EPPI::BCProps::calculate_primary_features () {
    // Set the size of the primary features matrix appropriately
    m_primaryFeatures.clear(); 
    m_primaryFeatures.resize(EPPI::PrimaryFeaturesList.size());
    // Set the base features
    for(size_t i=0; i<m_baseFeatures.size(); ++i) {
        m_primaryFeatures[i] = m_baseFeatures[i];}
    // Calculate the other primary features as appropriate
    m_primaryFeatures[21] = m_baseFeatures[0] + m_baseFeatures[1];
    m_primaryFeatures[22] = m_primaryFeatures[21] + m_baseFeatures[2];
    m_primaryFeatures[23] = m_baseFeatures[3] + m_baseFeatures[4];
    m_primaryFeatures[24] = m_primaryFeatures[23] + m_baseFeatures[5];
    m_primaryFeatures[25] = m_baseFeatures[6] + m_baseFeatures[7];
    m_primaryFeatures[26] = m_primaryFeatures[25] + m_baseFeatures[8];
    m_primaryFeatures[27] = m_baseFeatures[9] + m_baseFeatures[10];
    m_primaryFeatures[28] = m_primaryFeatures[27] + m_baseFeatures[11];
    m_primaryFeatures[29] = m_baseFeatures[15] + m_baseFeatures[16];
    m_primaryFeatures[30] = m_primaryFeatures[29] + m_baseFeatures[17];
    m_primaryFeatures[31] = m_baseFeatures[12] - m_primaryFeatures[30];
    m_primaryFeatures[32] = m_baseFeatures[18] + m_baseFeatures[19];
    m_primaryFeatures[33] = m_primaryFeatures[32] + m_baseFeatures[20];
    m_primaryFeatures[34] = m_primaryFeatures[22] + m_primaryFeatures[24];
    m_primaryFeatures[35] = m_primaryFeatures[26] + m_primaryFeatures[28];
    // Do an error checking to make sure the values are appropriate and
    // reasonable for an analysis
    error_check_primary ();
}

// Calculate and store a derived feature
void EPPI::BCProps::calculate_derived_feature (const Feature& which) {
    // Store the value here
    double value = 0.0;
    // Validate that the feature has content
    if (which.primary_features() == 0) {
        string error = "Error in calculate_derived_feature\n";
        throw PANTZ_error (error);}
    // Get the index of the first feature
    size_t I = primary_feature_index (which.primary_feature(0));
    value = m_primaryFeatures[I];
    // A flag for if divide by zero has happened
    bool divideFlag = false;
    // If there are further features
    if (which.primary_features() > 1) {
        for(size_t i=1; i<which.primary_features(); ++i) {
            // Get the index of the feature
            I = primary_feature_index(which.primary_feature(i));
            // Get the mathematical operation
            char L = which.operation(i-1);
            // Do the appropriate calculation
            if (L == '*') {value = value * m_primaryFeatures[I];}
            // Division needs special handling to take care of divide by 0 cases
            else if (L == '/') {
                // Get the magnitude of the primary feature
                double mag = fabs(m_primaryFeatures[I]);
                // If it is effectively 0
                if (mag < 1e-12) {
                    // If the numerator is also effectively 0, allow it
                    if (fabs(value) < 1e-12) {
                        // Just set the value to 0, and indicate that divide by
                        // 0 occurred
                        value = 0.0; 
                        divideFlag = true;}
                    // If the numerator is not 0, throw an error
                    else {
                        string error = "Divide by 0 error happened when "
                                       "calculating " + which.name() + " for "
                                     + m_name + "\n";
                        throw PANTZ_error (error);}}
                // If it is not effectively 0, do the calculation
                else {value = value / m_primaryFeatures[I];}}
            // Addition
            else if (L == '+') {value = value + m_primaryFeatures[I];}
            // Subtraction
            else if (L == '-') {value = value - m_primaryFeatures[I];}
            // Throw an error for any other operation
            else {
                string error = "Algorithm Error: calculate_derived_features "
                               "can't handle ";
                error += L;
                error += " as an operation.\n";
                throw PANTZ_error (error);}}}
    // If divide by zero happened, leave the value as 0
    if (divideFlag) {value = 0.0;}
    // Make sure the index of the feature matches where the value will be stored
    if (which.index() != m_derivedFeatures.size()) {
        string error = "Algorithm Error: calculating derived features out of "
                       "order.\n";
        throw PANTZ_error (error);}
    m_derivedFeatures.push_back(value);
}

// Calculate all of the derived features for a BCProp
void EPPI::BCProps::calculate_all_derived_features (const vector<Feature>& props) {
    // Clear existing features
    m_derivedFeatures.clear();
    // Allocate an appropriate amount of memory
    m_derivedFeatures.reserve(props.size());
    // Store each individual feature
    for(size_t i=0; i<props.size(); ++i) {
        calculate_derived_feature(props[i]);}
}
