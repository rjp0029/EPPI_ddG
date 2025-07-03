/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements methods of the BCProps class for finding the indexes of
 * base and primary features */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error indexes.h must be included by BCProps.h
#endif

// Implement the method for finding base indexes
size_t EPPI::BCProps::base_feature_index (const string& label) {
    // Loop through the base features
    for(size_t i=0; i<EPPI::BaseFeaturesList.size(); ++i) {
        if (label == EPPI::BaseFeaturesList[i]) {return i;}}
    // If the function reached this point, it is an unknown base feature
    string error = label + " is an unknown EPPI base feature.\n";
    throw PANTZ_error (error);
    // Return an invalid index just in case
    return EPPI::BaseFeaturesList.size();
}

// Same logic, but for primary features
size_t EPPI::BCProps::primary_feature_index (const string& label) {
    for(size_t i=0; i<EPPI::PrimaryFeaturesList.size(); ++i) {
        if (label == EPPI::PrimaryFeaturesList[i]) {return i;}}
    string error = label + " is an unknown EPPI primary feature.\n";
    throw PANTZ_error (error);
    return EPPI::PrimaryFeaturesList.size();
}
