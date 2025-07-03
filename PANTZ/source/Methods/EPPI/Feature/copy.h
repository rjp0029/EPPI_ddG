/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the copy method of the EPPI::Feature class */

// Make sure it is being loaded from the expected file
#ifndef EPPI_Feature_Loading_Status
#error EPPI::Feature::copy must be included by Feature.h
#endif

// Implement the function
void EPPI::Feature::copy (const Feature& other) {
    // Copy single values
    m_name = other.m_name;
    m_index = other.m_index;
    // Delete the primary vector's contents and resize it appropriately
    m_primary.clear(); m_primary.resize(other.m_primary.size());
    // If there is content, copy each value
    if (other.m_primary.size() > 0) {
        for(size_t i=0; i<other.m_primary.size(); ++i) {
            m_primary[i] = other.m_primary[i];}}
    // Use the same logic on the operations
    m_operations.clear(); m_operations.resize(other.m_operations.size());
    if (other.m_operations.size() > 0) {
        for(size_t i=0; i<other.m_operations.size(); ++i) {
            m_operations[i] = other.m_operations[i];}}
}
