/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the initialization method of the EPPI::Feature class */

// Make sure it is being loaded from the expected file
#ifndef EPPI_Feature_Loading_Status
#error EPPI::Feature::initialize must be included by Feature.h
#endif

// Implement the function
void EPPI::Feature::initialize () {
    m_name = "NULL";
    m_index = 0;
    m_primary.clear();
    m_operations.clear();
}
