/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the default initialization method of the NeighborBox
 * class */

// Make sure it is being included as expected
#ifndef NeighborBox_Loading_Status
#error The NeighborBox initialization file must be included by NeighborBox.h
#endif

// Implement the method
void PROT::NeighborBox::initialize () {
    // Assign default variable values
    m_XL = 0.0;
    m_XU = 0.0;
    m_YL = 0.0;
    m_YU = 0.0;
    m_ZL = 0.0;
    m_ZU = 0.0;
    m_bin_side = 10000.0;
    m_Xbins = 0;
    m_Ybins = 0;
    m_Zbins = 0;
    m_size = 0;
}
