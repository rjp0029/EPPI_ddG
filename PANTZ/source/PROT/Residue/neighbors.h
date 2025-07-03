/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * rotamer methods for a residue
*/

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

#include "../Protein.h"
#include "../KDtree.h"

// get the intra protein neighbors of this residue
vector<PROT::Residue*> PROT::Residue::get_intra_neighbors(PROT::Protein* protein, float cutoff) {
    // get the residues in the protein
    vector<PROT::Residue*> residues;
    for (size_t i = 0; i < protein->size(); i++) {
        residues.push_back(protein->operator()(i, ' ', true));
    }
    // get the kdtree of the protein
    KDtree<PROT::Residue> tree(residues);
    // get the neighbors of this residue
    // m_intra_neighbors = tree.radius_neighbors(this, cutoff);
    return tree.radius_neighbors(this, cutoff);
}

// get the interprotein neighbors of this residue
vector<PROT::Residue*> PROT::Residue::get_inter_neighbors(vector<PROT::Protein*> proteins, float cutoff) {
    // get the residues in the protein
    vector<PROT::Residue*> residues;
    for (size_t i = 0; i < proteins.size(); i++) {
        for (size_t j = 0; j < proteins[i]->size(); j++) {
            residues.push_back(proteins[i]->operator()(j, ' ', true));
        }
    }
    // get the kdtree of the protein
    KDtree<PROT::Residue> tree(residues);
    // get the neighbors of this residue
    // m_inter_neighbors = tree.radius_neighbors(this, cutoff);
    return tree.radius_neighbors(this, cutoff);
}

// get neighbors from a vector of residues 
vector<PROT::Residue*> PROT::Residue::get_inter_neighbors_res(vector<PROT::Residue*> residues, float cutoff) {
    // get the kdtree of the protein
    KDtree<PROT::Residue> tree(residues);
    // get the neighbors of this residue
    // m_neighbors = tree.radius_neighbors(this, cutoff);
    return tree.radius_neighbors(this, cutoff);
}

// get intra neighbors from a vector of residues
vector<PROT::Residue*> PROT::Residue::get_intra_neighbors(vector<PROT::Residue*> residues, float cutoff) {
    // get the kdtree of the protein
    KDtree<PROT::Residue> tree(residues);
    // get the neighbors of this residue
    // m_neighbors = tree.radius_neighbors(this, cutoff);
    return tree.radius_neighbors(this, cutoff);
}