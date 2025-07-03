/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the calculate_bin method of the NeighborBox class. This
 * function calculates the X, Y, or Z bin of a specified point */

// Make sure it is being included as expected
#ifndef NeighborBox_Loading_Status
#error The NeighborBox calculate file must be included by NeighborBox.h
#endif

// Implement the method
size_t PROT::NeighborBox::calculate_bin (const coor where, const char L) {
    // The answer from the function
    size_t answer = 0;
    // Identify the appropriate lower bound to use for the calculation
    coor LB;
    if (L == 'X') {LB = m_XL;}
    else if (L == 'Y') {LB = m_YL;}
    else {LB = m_ZL;}
    // If the coordinate is less than the lower bound, return 0
    if (where < LB) {answer = 0; return answer;}
    // Calculate the bin index
    answer = (size_t) ((where - LB) / m_bin_side);
    // If it exceeds the bounds, reset it to the maximum possible value
    if (L == 'X') {if (answer >= m_Xbins) {answer = m_Xbins-1;}}
    else if (L == 'Y') {if (answer >= m_Ybins) {answer = m_Ybins-1;}}
    else {if (answer >= m_Zbins) {answer = m_Zbins-1;}}
    // Return the answer
    return answer;
}
