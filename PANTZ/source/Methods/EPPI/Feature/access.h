/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains methods of the EPPI::Feature class that allow access to
 * its class attributes.*/

// Make sure it is being loaded from the expected file
#ifndef EPPI_Feature_Loading_Status
#error EPPI::Feature::access.h must be included by Feature.h
#endif

// Access to a primary feature
string EPPI::Feature::primary_feature (const size_t i) const {
    // Make sure i is a valid index
    if (i >= m_primary.size()) {
        stringstream c1; c1 << m_primary.size();
        stringstream c2; c2 << i;
        string error = m_name + " contains " + c1.str() + " primary features, so "
                     + c2.str() + " is an invalid index\n";
        throw PANTZ_error (error);}
    // Return the value
    return m_primary[i];
}

// Access to a math operation
char EPPI::Feature::operation (const size_t i) const {
    if (i >= m_operations.size()) {
        stringstream c1; c1 << m_operations.size();
        stringstream c2; c2 << i;
        string error = m_name + " contains " + c1.str() + " operations, so "
                     + c2.str() + " is an invalid index\n";
        throw PANTZ_error (error);}
    return m_operations[i];
}
