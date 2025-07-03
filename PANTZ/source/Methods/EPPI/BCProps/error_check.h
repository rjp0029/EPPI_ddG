/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the method of the BCProps class responsible for error
 * checking the base and derived features to make sure they're realistic */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error error_check.h must be included by BCProps.h
#endif

void EPPI::BCProps::error_check_primary () {
    // If the primary features list isn't the right size, throw an error
    if (m_primaryFeatures.size() != PrimaryFeaturesList.size()) {
        string error = "Algorithm Error: error checking a BCProps with an "
                       "improper number of primary features\n";
        throw PANTZ_error (error);}
    // Store a list of primary features that have issues here
    vector<string> issues;
    // Check each feature. This is hard coded :(
    string BOUNDS = "-0.5 1000"  // n_epp_sb
        " -0.5 1000"             // n_epp_hb
        " -0.5 1000"             // n_epp_hp
        " -0.5 1000"             // n_np_sb
        " -0.5 1000"             // n_np_hb
        " -0.5 1000"             // n_np_hp
        " -0.5 1000000"          // rl_epp_sb
        " -0.5 1000000"          // rl_epp_hb
        " -0.5 1000000"          // rl_epp_hp
        " -0.5 1000000"          // rl_np_sb
        " -0.5 1000000"          // rl_np_hb
        " -0.5 1000000"          // rl_np_hp
        " -10000 10000"           // dg_separated
        " 50 1000000"            // bsa
        " 0 1"                   // sc
        " -10000 10000"           // epbe_sb
        " -10000 10000"           // epbe_hb
        " -10000 10000"           // epbe_hp
        " -0.5 1000"             // n_st_st_sb
        " -0.5 1000"             // n_st_st_hb
        " -0.5 1000"             // n_st_st_hp
        " -0.5 2000"             // n_epp_polar
        " 0.5 3000"              // n_eppi
        " -0.5 2000"             // n_np_polar
        " -0.5 3000"             // n_nppi
        " -0.5 2000000"          // rl_epp_polar
        " 0.5 3000000"           // rl_eppi
        " -0.5 2000000"          // rl_np_polar
        " -0.5 3000000"          // rl_nppi
        " -10000 2000"           // epbe_polar
        " -10000 3000"           // epbe
        " -10000 3000"           // npbe
        " -0.5 2000"             // n_st_st_polar
        " -0.5 3000"             // n_st_st_pi
        " 0.5 5000"              // n_total
        " 0.5 5000000";          // rl_total
    // Split the bounds into a vector of strings
    vector<string> bounds; Text::split(bounds, BOUNDS);
    // A list of which features are problematic
    vector<string> problems;
    // Do this for each property
    for (size_t i=0; i<PrimaryFeaturesList.size(); ++i) {
        // The two indices to check
        size_t i1 = 2*i;
        size_t i2 = i1 + 1;
        // Convert those values to numbers
        double b1; stringstream c1; c1 << bounds[i1]; c1 >> b1;
        double b2; stringstream c2; c2 << bounds[i2]; c2 >> b2;
        // Compare them to the value
        if ((m_primaryFeatures[i] < b1) || (m_primaryFeatures[i] > b2)) {
            problems.push_back(PrimaryFeaturesList[i]);}}
    // If there are primary features that are a problem, throw an error
    if (problems.size() > 0) {
        string error = "Binding Complex: " + m_name + "\nProblematic Feature";
        if (problems.size() > 1) {error += "s";}
        error += ":";
        for(size_t i=0; i<problems.size(); ++i) {
            // Add the feature
            error += " " + problems[i];
            // Possibly add a comma
            if ((problems.size() > 2) && (i < problems.size()-1)) {error += ",";}
            // Possibly add "and"
            if ((problems.size() > 1) && (i == problems.size()-2)) {
                error += " and";}}
        error += "\n";
        // List out the primary features
        for(size_t i=0; i<PrimaryFeaturesList.size(); ++i) {
            stringstream c1; c1 << fixed << setprecision(4) << m_primaryFeatures[i];
            string label = PrimaryFeaturesList[i] + ":";
            Text::ljust_insert(error, label, 15, ' ');
            error += c1.str() + "\n";}
        // Throw the error
        throw PANTZ_error (error);}
}
