/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to find hydrogen bonds between two residues. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

#include <map>

map<string, vector<string>> polar_atoms_charmm = {
        {"GLY", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"ALA", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"VAL", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"CYS", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"PRO", {"O", "C", "OT1", "OT2", "HT1", "HT2", "HT3", "N"}},
        {"LEU", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"ILE", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"MET", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"TRP", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "HE1", "N", "NE1"}},
        {"PHE", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "N"}},
        {"SER", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OG", "HG", "N"}},
        {"THR", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OG1", "HG1", "N"}},
        {"TYR", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OH", "HH", "N"}},
        {"ASN", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OD1", "HD21", "HD22", "N", "ND2"}},
        {"GLN", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OE1", "HE21", "HE22", "N", "NE2"}},
        {"LYS", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "HZ1", "HZ2", "HZ3", "N", "NZ"}},
        {"ARG", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "HE", "HH11", "HH12", "HH21", "HH22", "N", "NH1", "NH2", "NE"}},
        {"HIS", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "ND1", "HE2", "NE2", "HD1", "N"}},
        {"HSD", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "ND1", "HE2", "NE2", "HD1", "N"}},
        {"ASP", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OD1", "OD2", "N"}},
        {"GLU", {"O", "C", "OT1", "OT2", "HN", "HT1", "HT2", "HT3", "OE1", "OE2", "N"}}
    };

// bool to determine if the passed atom is polar in this residue
bool is_polar(const PROT::Atom* atom) {
    // get the residue name
    string res_name = atom->residue();
    // get the atom name
    string atom_name = atom->name();
    // strip the atom name
    Text::strip(atom_name);
    // find the name in polar atoms charmm map map<string, vector<string>> polar_atoms_charmm
    // if the atom is in the polar atoms map
    if (find(polar_atoms_charmm.at(res_name).begin(), polar_atoms_charmm.at(res_name).end(), atom_name) != polar_atoms_charmm.at(res_name).end()) {
        // return true
        return true;
    }
    // }
    // otherwise return false
    return false;
}

// determine the number of hydrophobic interactions between two residues
vector<PROT::Hydrophobic> PROT::Residue::hydrophobic(PROT::Residue* other, float sasa_cutoff, bool verbose) {
    vector<PROT::Hydrophobic> hydrophobic_interactions;
    // this assumes that the phobic sasa points have been set
    // the water probe radius
    float probe = 1.4;
    // number of mesh points, Shrake and Rupley used 92, Varun used 80
    size_t n_mesh = 80;
    // the buried surface area of the residue to residue interaction with all combinations of backbone and side chain interactions
    float bb_bb_bsa_1 = 0;
    float bb_bb_bsa_2 = 0;
    float bb_sc_bsa_1 = 0;
    float bb_sc_bsa_2 = 0;
    float sc_bb_bsa_1 = 0;
    float sc_bb_bsa_2 = 0;
    float sc_sc_bsa_1 = 0;
    float sc_sc_bsa_2 = 0;

    // go through all residues and get the sc_sc_bsa
    // for (size_t i = 0; i < m_count; i++) {
        // PROT::Atom* atom = &m_atoms[i];
    for (size_t i = 0; i < m_count; i++) {
        PROT::Atom* atom = &m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            // for (size_t j = 0; j < other->m_count; j++) {
                // PROT::Atom* atom2 = &other->m_atoms[j];
            for (size_t j = 0; j < other->m_count; j++) {
                PROT::Atom* atom2 = &other->m_atoms[j];
                // if the atom is a backbone atom, skip
                if (atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        sc_sc_bsa_1 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }
    
    // go through all atoms and get the sc_sc_bsa 2
    // for (size_t i = 0; i < other->m_count; i++) {
        // PROT::Atom* atom = other->&m_atoms[i];
    for (size_t i = 0; i < other->m_count; i++) {
        PROT::Atom* atom = &other->m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            // for (size_t j = 0; j < m_count; j++) {
                // PROT::Atom* atom2 = &m_atoms[j];
            for (size_t j = 0; j < m_count; j++) {
                PROT::Atom* atom2 = &m_atoms[j];
                // if bb skip
                if (atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    // cout<<"bp: "<<other->m_name<<' '<<atom.m_name<<' '<<point[0]<<' '<<point[1]<<' '<<point[2]<<" by "<<m_name<<' '<<atom2.m_name<<' '<<bsa2<<endl;
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        sc_sc_bsa_2 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the bb_sc_bsa 1
    for (size_t i = 0; i < m_count; i++) {
        PROT::Atom* atom = &m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (!atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < other->m_count; j++) {
                PROT::Atom* atom2 = &other->m_atoms[j];
                // if bb skip
                if (atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        bb_sc_bsa_1 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the bb_sc_bsa 2
    for (size_t i = 0; i < other->m_count; i++) {
        PROT::Atom* atom = &other->m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < m_count; j++) {
                PROT::Atom* atom2 = &m_atoms[j];
                // if bb skip
                if (!atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        bb_sc_bsa_2 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the sc_bb_bsa 1
    for (size_t i = 0; i < m_count; i++) {
        PROT::Atom* atom = &m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < other->m_count; j++) {
                PROT::Atom* atom2 = &other->m_atoms[j];
                // if bb skip
                if (!atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        sc_bb_bsa_1 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the sc_bb_bsa 2
    for (size_t i = 0; i < other->m_count; i++) {
        PROT::Atom* atom = &other->m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (!atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < m_count; j++) {
                PROT::Atom* atom2 = &m_atoms[j];
                // if bb skip
                if (atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        sc_bb_bsa_2 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the bb_bb_bsa 1
    for (size_t i = 0; i < m_count; i++) {
        PROT::Atom* atom = &m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (!atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < other->m_count; j++) {
                PROT::Atom* atom2 = &other->m_atoms[j];
                // if bb skip
                if (!atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        bb_bb_bsa_1 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through all atoms and get the bb_bb_bsa 2
    for (size_t i = 0; i < other->m_count; i++) {
        PROT::Atom* atom = &other->m_atoms[i];
        // if the atom is polar, skip
        if (is_polar(atom)) {
            continue;
        }
        // if backbone atom, skip
        if (!atom->is_backbone_atom()) {
            continue;
        }
        // go through the atom's sa points and determine if they are buried
        size_t n_solv = 0;
        for (auto & point : atom->m_sasa_points) {
            bool solvent_accessible = true;
            for (size_t j = 0; j < m_count; j++) {
                PROT::Atom* atom2 = &m_atoms[j];
                // if bb skip
                if (!atom2->is_backbone_atom()) {
                    continue;
                }
                // get the distance between the mesh point and atom2 
                float dist = sqrt(pow(point[0] - atom2->m_coors[0], 2) + 
                                  pow(point[1] - atom2->m_coors[1], 2) + 
                                  pow(point[2] - atom2->m_coors[2], 2));
                // if the distance is less than the radius of atom2
                if (dist < atom2->lj_sigma() + probe) {
                    solvent_accessible = false;
                    break;
                }
            }
            if (!solvent_accessible) {
                n_solv++;
            }
        }
        // calculate the effective radius of the atom
        float effrad = atom->lj_sigma() + probe;
        // add the solvent accessible surface area of the atom to the bsa
        bb_bb_bsa_2 += ((float(n_solv) / float(n_mesh)) * 4 * M_PI * effrad * effrad);
    }

    // go through the different sums of bsa 1 and 2 and if any sums are greater 
    // than sasa_cutoff then there is a hydrophobic interaction, return the number 
    // of sums greater than sasa_cutoff and the type
    float bb_bb_bsa = bb_bb_bsa_1 + bb_bb_bsa_2;
    float bb_sc_bsa = bb_sc_bsa_1 + bb_sc_bsa_2;
    float sc_bb_bsa = sc_bb_bsa_1 + sc_bb_bsa_2;
    float sc_sc_bsa = sc_sc_bsa_1 + sc_sc_bsa_2;

    // get the number of bsas that are greater than sasa_cutoff
    size_t count = 0;
    if (bb_bb_bsa > sasa_cutoff) {
        PROT::Hydrophobic interaction;
        interaction.res1 = this;
        interaction.res1_backbone = true;
        interaction.res2 = other;
        interaction.res2_backbone = true;
        interaction.phobic_BSASA = bb_bb_bsa;
        hydrophobic_interactions.push_back(interaction);
    }
    if (bb_sc_bsa > sasa_cutoff) {
        PROT::Hydrophobic interaction;
        interaction.res1 = this;
        interaction.res1_backbone = true;
        interaction.res2 = other;
        interaction.res2_backbone = false;
        interaction.phobic_BSASA = bb_sc_bsa;
        hydrophobic_interactions.push_back(interaction);
    }
    if (sc_bb_bsa > sasa_cutoff) {
        PROT::Hydrophobic interaction;
        interaction.res1 = this;
        interaction.res1_backbone = false;
        interaction.res2 = other;
        interaction.res2_backbone = true;
        interaction.phobic_BSASA = sc_bb_bsa;
        hydrophobic_interactions.push_back(interaction);
    }
    if (sc_sc_bsa > sasa_cutoff) {
        PROT::Hydrophobic interaction;
        interaction.res1 = this;
        interaction.res1_backbone = false;
        interaction.res2 = other;
        interaction.res2_backbone = false;
        interaction.phobic_BSASA = sc_sc_bsa;
        hydrophobic_interactions.push_back(interaction);
    }
    if ((verbose) and (hydrophobic_interactions.size() > 0)) {
        cout << "Hydrophobic interactions found:\n";
        for (auto & interaction : hydrophobic_interactions) {
            cout << interaction.res1->m_name << ' ' << interaction.res1->m_number << ' ' << interaction.res1->m_protein << ' ';
            if (interaction.res1_backbone) {
                cout << "bb ";
            } else {
                cout << "sc ";
            }
            cout << interaction.res2->m_name << ' ' << interaction.res2->m_number << ' ' << interaction.res2->m_protein << ' ';
            if (interaction.res2_backbone) {
                cout << "bb ";
            } else {
                cout << "sc ";
            }
            cout << interaction.phobic_BSASA << " angstroms^2\n\n";
        }
    }
    return hydrophobic_interactions;
}
