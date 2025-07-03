/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the initialization method of the BCProps class */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error initialization.h must be included by BCProps.h
#endif

// Implement the method
void EPPI::BCProps::initialize () {
    m_name = "NULL";
    m_baseFeatures.clear();
    m_primaryFeatures.clear();
    m_derivedFeatures.clear();
}
