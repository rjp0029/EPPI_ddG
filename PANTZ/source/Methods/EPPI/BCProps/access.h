/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the operator[] method of the BCProps class */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error access.h must be included by BCProps.h
#endif

// Implement the method
double EPPI::BCProps::operator[] (const size_t i) const {
    // Validate the index
    if (i >= m_derivedFeatures.size()) {
        stringstream c1; c1 << m_derivedFeatures.size();
        stringstream c2; c2 << i;
        string error = "The " + m_name + " Binding Complex contains "
                     + c1.str() + " derived features, so " + c2.str()
                     + " is an invalid index.\n";
        throw PANTZ_error (error);}
    // Return the value
    return m_derivedFeatures[i];
}
