/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the copy method of the BCProps class */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error copy.h must be included by BCProps.h
#endif

// Implement the method
void EPPI::BCProps::copy (const BCProps& other) {
    // Copy the name
    m_name = other.m_name;
    // Clear the base features and reallocate memory
    m_baseFeatures.clear(); m_baseFeatures.reserve(other.m_baseFeatures.size());
    // If there are values, copy them
    if (other.m_baseFeatures.size() > 0) {
        for(size_t i=0; i<other.m_baseFeatures.size(); ++i) {
            m_baseFeatures.push_back(other.m_baseFeatures[i]);}}
    // follow the same logic for the primary features
    m_primaryFeatures.clear(); m_primaryFeatures.reserve(other.m_primaryFeatures.size());
    if (other.m_primaryFeatures.size() > 0) {
        for(size_t i=0; i<other.m_primaryFeatures.size(); ++i) {
            m_primaryFeatures.push_back(other.m_primaryFeatures[i]);}}
    // Follow the same logic for the derived features
    m_derivedFeatures.clear(); m_derivedFeatures.reserve(other.m_derivedFeatures.size());
    if (other.m_derivedFeatures.size() > 0) {
        for(size_t i=0; i<other.m_derivedFeatures.size(); ++i) {
            m_derivedFeatures.push_back(other.m_derivedFeatures[i]);}}
}
