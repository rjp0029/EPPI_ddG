/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to find a salt bridge between two residues. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif



// function to find a salt bridge with another residue
vector<PROT::SaltBridge> PROT::Residue::salt_bridge(PROT::Residue* other, float distance, bool verbose) {
    // combinations possible are ARG and ASP, ARG and GLU, ARG and c_term,
    // LYS and ASP, LYS and GLU, LYS and c_term, n_term and ASP, n_term and GLU
    // n_term and c_term
    if (!((m_name == "ARG" && (other->m_name == "ASP" || other->m_name == "GLU" || other->m_C_terminus)) ||
          (m_name == "LYS" && (other->m_name == "ASP" || other->m_name == "GLU" || other->m_C_terminus)) ||
          (m_N_terminus && (other->m_name == "ASP" || other->m_name == "GLU" || other->m_C_terminus)) ||
          (other->m_name == "ARG" && (m_name == "ASP" || m_name == "GLU" || m_C_terminus)) ||
          (other->m_name == "LYS" && (m_name == "ASP" || m_name == "GLU" || m_C_terminus)) ||
          (other->m_N_terminus && (m_name == "ASP" || m_name == "GLU" || m_C_terminus)))) {
        // if none of these are true then return empty map
        return {};
    }
    vector<PROT::SaltBridge> salt_bridges;
    // Loop through the atoms in the residue
    // for (PROT::Atom& atom : m_atoms) {
    for (size_t i_atom = 0; i_atom < m_count; i_atom++){
        PROT::Atom* atom = &m_atoms[i_atom];
        // Check if the atom is a nitrogen
        if (atom->m_element == "N") {
            // if this is a backbone nitrogen, and the residue is not a C-terminus, skip
            if (atom->is_backbone_atom() && !m_C_terminus) {
                continue;
            }
            // Loop through the atoms in the other residue
            // for (PROT::Atom& other_atom : other->m_atoms) {
            for (size_t j_atom = 0; j_atom < other->m_count; j_atom++){
                PROT::Atom* other_atom = &other->m_atoms[j_atom];
                // Check if the atom is an oxygen or sulfur
                if (other_atom->m_element == "O" || other_atom->m_element == "S") {
                    // if this is a backbone oxygen, and the residue is not a N-terminus, skip
                    if (other_atom->is_backbone_atom() && !other->m_N_terminus) {
                        continue;
                    }
                    // Calculate squared distance between atoms
                    float dist = atom->distance(*other_atom);
                    if (dist < distance) {
                        if (verbose) {
                            cout << "Salt bridge found with distance: "<<dist<<"\n";
                            cout << atom->str();
                            cout << other_atom->str()<<"\n";
                        }
                        PROT::SaltBridge salt_bridge;
                        salt_bridge.donor_residue = this;
                        salt_bridge.donor_backbone = atom->is_backbone_atom();
                        salt_bridge.acceptor_residue = other;
                        salt_bridge.acceptor_backbone = other_atom->is_backbone_atom();
                        salt_bridge.donor_atom = atom;
                        salt_bridge.acceptor_atom = other_atom;
                        salt_bridge.distance = dist;
                        salt_bridges.push_back(salt_bridge);
                    }
                }
            }
        }
    }
    // Loop through the atoms in the other residue to check for the reverse salt bridge
    // for (PROT::Atom& other_atom : other->m_atoms) {
    for (size_t j_atom = 0; j_atom < other->m_count; j_atom++){
        PROT::Atom * other_atom = &other->m_atoms[j_atom];
        // Check if the atom is a nitrogen
        if (other_atom->m_element == "N") {
            // if this is a backbone nitrogen, and the residue is not a C-terminus, skip
            if (other_atom->is_backbone_atom() && !other->m_C_terminus) {
                continue;
            }
            // Loop through the atoms in the other residue
            // for (PROT::Atom& atom : m_atoms) {
            for (size_t i_atom = 0; i_atom < m_count; i_atom++){
                PROT::Atom * atom = &m_atoms[i_atom];
                // Check if the atom is an oxygen or sulfur
                if (atom->m_element == "O" || atom->m_element == "S") {
                    // if this is a backbone oxygen, and the residue is not a N-terminus, skip
                    if (atom->is_backbone_atom() && !m_N_terminus) {
                        continue;
                    }
                    // Calculate squared distance between atoms
                    float dist = other_atom->distance(*atom);
                    if (dist < distance) {
                        if (verbose) {
                            cout << "Salt bridge found with distance: "<<dist<<"\n";
                            cout << atom->str();
                            cout << other_atom->str()<<"\n";
                        }
                        PROT::SaltBridge salt_bridge;
                        salt_bridge.donor_residue = other;
                        salt_bridge.donor_backbone = other_atom->is_backbone_atom();
                        salt_bridge.acceptor_residue = this;
                        salt_bridge.acceptor_backbone = atom->is_backbone_atom();
                        salt_bridge.donor_atom = other_atom;
                        salt_bridge.acceptor_atom = atom;
                        salt_bridge.distance = dist;
                        salt_bridges.push_back(salt_bridge);
                    }
                }
            }
        }
    }
    // remove duplicate salt bridges that are between the same residues and backbone types
    // vector to store the unique salt bridges
    vector<PROT::SaltBridge> unique_salt_bridges;
    // loop through the salt bridges
    for (PROT::SaltBridge& salt_bridge : salt_bridges) {
        // if the salt bridge is not in the unique salt bridges
        if (find(unique_salt_bridges.begin(), unique_salt_bridges.end(), salt_bridge) == unique_salt_bridges.end()) {
            // add the salt bridge to the unique salt bridges
            unique_salt_bridges.push_back(salt_bridge);
        }
    }
    return unique_salt_bridges;
}
