/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements the private copy method of the class, which is used for
 * copying information from one instance of the class to another. */ 

// Make sure it is being included as expected
#ifndef NeighborBox_Loading_Status
#error The NeighborBox copy file must be included by NeighborBox.h
#endif

// Implement the method
void PROT::NeighborBox::copy (const NeighborBox& other) {
    // Copy the non-vector values
    m_XL = other.m_XL;
    m_XU = other.m_XU;
    m_YL = other.m_YL;
    m_YU = other.m_YU;
    m_ZL = other.m_ZL;
    m_ZU = other.m_ZU;
    m_bin_side = other.m_bin_side;
    m_Xbins = other.m_Xbins;
    m_Ybins = other.m_Ybins;
    m_Zbins = other.m_Zbins;
    m_size = other.m_size;
    // Copy the vector of atoms
    m_atom_pointers.clear();
    m_atom_pointers.resize(m_Xbins);
    for(size_t i=0; i<m_Xbins; ++i) {
        m_atom_pointers[i].resize(m_Ybins);
        for(size_t j=0; j<m_Ybins; ++j) {
            m_atom_pointers[i][j].resize(m_Zbins);
            for(size_t k=0; k<m_Zbins; ++k) {
                m_atom_pointers[i][j][k].resize(other.m_atom_pointers[i][j][k].size());
                for(size_t l=0; l<other.m_atom_pointers[i][j][k].size(); ++l) {
                    m_atom_pointers[i][j][k][l] = other.m_atom_pointers[i][j][k][l];
                }}}}
    // End the function
}
