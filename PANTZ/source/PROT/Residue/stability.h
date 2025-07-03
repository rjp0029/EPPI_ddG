/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h 
    file. It includes methods associated with setting the stbaility 
    of a rotamer residue both in the context of prestabilitzation
    (one protein) and bound stabilization (vector of proteins).
*/

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

#include "../Protein.h"

size_t PROT::Residue::free_rotamers(vector<PROT::Residue> rotamers, vector<PROT::Residue*> neighbors) {
    // count the clashes with the neighbors
    size_t count = 0;
    for (size_t i = 0; i < rotamers.size(); i++) {
        for (size_t j = 0; j < neighbors.size(); j++) {
            if (rotamers[i].heavy_side_chain_clash(neighbors[j])) {
            // if (rotamers[i].heavy_clash(neighbors[j])) {
                count++;
                // skip to the next rotamer
                break;
            }
        } 
    }
    // return the number of rotamers that are free of clashes
    return rotamers.size() - count;
}

// function to check the stability of the residue in its environment
void PROT::Residue::set_stability(vector<PROT::Residue> rotamers, vector<PROT::Residue*> neighbors, bool backbone, string& how) {
    // if backbone set to true
    if (backbone) {
        if (how == "pre") {
            m_prestable = true;
        }
        else {
            m_bound_stable = true;
        }
        return;
    }
    if (m_name == "GLY" || m_name == "ALA" || m_name == "PRO") {
        if (how == "pre") {
            m_prestable = true;
        }
        else {
            m_bound_stable = true;
        }
        return;
    }
    bool salt_bridge_possible = false;
    size_t sb_count = 0;
    // check if a salt bridge is possible (name of residue is ARG or LYS or GLU or ASP or m_C_terminal or m_N_terminal)
    if (m_name == "ARG" || m_name == "LYS" ||
        m_name == "GLU" || m_name == "ASP" || 
        m_C_terminus == true || m_N_terminus == true) {
        salt_bridge_possible = true;
    }
    // check for salt bridge with neighbors
    if (salt_bridge_possible) {
        vector<PROT::SaltBridge> salt_bridges;
        for (size_t i = 0; i < neighbors.size(); i++) {
            salt_bridges = this->salt_bridge(neighbors[i], 4.0, false);
            if (salt_bridges.size() > 0) {
                sb_count += salt_bridges.size();
            }
        }
        // if sb > 1 set to prestable and stop checking
        if (sb_count > 1) {
            if (how == "pre") {
                m_prestable = true;
            }
            else {
                m_bound_stable = true;
            }
            return;
        }
    }
    // get the number of free rotamers
    size_t free_rot = free_rotamers(rotamers, neighbors);
    // check if prestable
    if ((rotamers.size() - free_rot) > 0.8*rotamers.size()) {
        if (how == "pre") {
            m_prestable = true;
        }
        else {
            m_bound_stable = true;
        }
    }
    else if ((m_name == "SER" || m_name == "THR" || m_name == "CYS" || m_name == "VAL") and free_rot < 3) {
        if (how == "pre") {
            m_prestable = true;
        }
        else {
            m_bound_stable = true;
        }
    }
    else {
        if (how == "pre") {
            m_prestable = false;
        }
        else {
            m_bound_stable = false;
        }
    }
}
