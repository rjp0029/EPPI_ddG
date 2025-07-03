/* Created by clay at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * method to find hydrogen bonds between two residues. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// Include the necessary header files
#include <map>
#include "../Atom.h"

// store the hbond donor and acceptor atoms for each residue: (From CHARMM topology files top_all36_prot.rtf)
map<char, vector<string>> hbond_donors = {
    {'A', {"HN N"}},
    {'R', {"HN N", "HE NE", "HH11 NH1", "HH12 NH1", "HH21 NH2", "HH22 NH2"}},
    {'N', {"HN N", "HD21 ND2", "HD22 ND2"}},
    {'D', {"HN N"}},
    {'C', {"HN N", "HG1 SG"}},
    {'Q', {"HN N", "HE21 NE2", "HE22 NE2"}},
    {'E', {"HN N"}},
    {'G', {"HN N"}},
    {'H', {"HN N", "HD1 ND1", "HE2 NE2"}},
    {'I', {"HN N"}},
    {'L', {"HN N"}},
    {'K', {"HN N", "HZ1 NZ", "HZ2 NZ", "HZ3 NZ"}},
    {'M', {"HN N"}},
    {'F', {"HN N"}},
    {'P', {}},
    {'S', {"HN N", "HG1 OG"}},
    {'T', {"HN N", "HG1 OG1"}},
    {'W', {"HN N", "HE1 NE1"}},
    {'Y', {"HN N", "HH OH"}},
    {'V', {"HN N"}}
};

map<char, vector<string>> hbond_acceptors = {
    {'A', {"O C"}},
    {'R', {"O C"}},
    {'N', {"OD1 CG", "O C"}},
    {'D', {"OD1 CG", "OD2 CG", "O C"}},
    {'C', {"O C"}},
    {'Q', {"OE1 CD", "O C"}},
    {'E', {"OE1 CD", "OE2 CD", "O C"}},
    {'G', {"O C"}},
    {'H', {"NE2 CD2_CE1", "O C", "ND1 CG_CE1"}},
    {'I', {"O C"}},
    {'L', {"O C"}},
    {'K', {"O C"}},
    {'M', {"O C"}},
    {'F', {"O C"}},
    {'P', {"O C"}},
    {'S', {"OG CB", "O C"}},
    {'T', {"OG1 CB", "O C"}},
    {'W', {"O C"}},
    {'Y', {"OH CZ", "O C"}},
    {'V', {"O C"}}
};

// function to count the number of hydrogen bonds between two residues
vector<PROT::HydrogenBond> PROT::Residue::hbond(Residue * other, float distance, float dha_angle, float daa_angle, bool verbose) {
    // initialize the count of hbonds to 0
    vector<HydrogenBond> hbonds;
    // get this residue's 1 letter code
    char this_residue = AA1();
    // get the other residue's 1 letter code
    char other_residue = other->AA1();
    // go through the atoms of this residue if they are hydrogens, get their donor atoms and then search for acceptor atoms
    for (size_t i_atom = 0; i_atom < m_count; i_atom++){
    // for (PROT::Atom hydrogen : m_atoms) {
        PROT::Atom* hydrogen = &m_atoms[i_atom];
        if (hydrogen->m_element != "H") {
            continue;
        }
        // get the donor strings 
        vector<string> donors = hbond_donors[this_residue];
        // go through donors and if the atom is a donatable hydrogen, get the donor atom
        for (string donor_pair : donors) {
            // split the donor pair, the donor is the second element
            vector<string> hydrogen_and_donor= Text::split(donor_pair, ' ');
            if (Text::contains(hydrogen->name(), hydrogen_and_donor[0])) {
                string donor_name = hydrogen_and_donor[1];
                Text::strip(donor_name);
                // get the donor atom by reference
                PROT::Atom* donor = this->get_atom(donor_name);
                // go through the atoms of the other residue and if they are acceptors, get their acceptor atoms
                for (size_t j_atom = 0; j_atom < other->m_count; j_atom++){
                    PROT::Atom* acceptor = &other->m_atoms[j_atom];
                    if (acceptor->m_element != "O" && acceptor->m_element != "N") {
                        continue;
                    }
                    bool two_antecedents = false;
                    // get the acceptor strings
                    vector<string> acceptors = hbond_acceptors[other_residue];
                    // go through acceptors and if the atom is an acceptor, get the acceptor atom
                    for (string acceptor_pair : acceptors) {
                        // split the acceptor pair, the acceptor is the first element
                        vector<string> acceptor_and_antecedent = Text::split(acceptor_pair, ' ');

                        // if the second element has an underscore, then there are 2 antecedents, get their midpoint as the antecedent
                        if (Text::contains(acceptor_and_antecedent[1], "_")) {
                            two_antecedents = true;
                        }
                        if (Text::contains(acceptor->name(), acceptor_and_antecedent[0])) {
                            PROT::Atom dummy;
                            PROT::Atom* antecedent;
                            // if there are two antecedents, get the midpoint
                            if (two_antecedents) {
                                string antecedent_name = acceptor_and_antecedent[1];
                                // split the antecedent pair
                                vector<string> antecedents = Text::split(antecedent_name, '_');
                                // get the antecedent atoms by reference
                                Text::strip(antecedents[0]);
                                Text::strip(antecedents[1]);
                                PROT::Atom * antecedent1 = other->get_atom(antecedents[0]);
                                PROT::Atom * antecedent2 = other->get_atom(antecedents[1]);
                                // get the midpoint of the two antecedents
                                dummy.m_coors[0] = (antecedent1->m_coors[0] + antecedent2->m_coors[0]) / 2;
                                dummy.m_coors[1] = (antecedent1->m_coors[1] + antecedent2->m_coors[1]) / 2;
                                dummy.m_coors[2] = (antecedent1->m_coors[2] + antecedent2->m_coors[2]) / 2;
                                dummy.m_name = "dummy";
                                antecedent = &dummy;
                            } else {
                                string antecedent_name = acceptor_and_antecedent[1];
                                Text::strip(antecedent_name);
                                antecedent = other->get_atom(antecedent_name);
                            }
                            // get the angle from the donor - hydrogen - acceptor using the dot product
                            float DHx = donor->x() - hydrogen->x();
                            float DHy = donor->y() - hydrogen->y();
                            float DHz = donor->z() - hydrogen->z();
                            float HAx = acceptor->x() - hydrogen->x();
                            float HAy = acceptor->y() - hydrogen->y();
                            float HAz = acceptor->z() - hydrogen->z();
                            float dot = DHx * HAx + DHy * HAy + DHz * HAz;
                            float DHmag = sqrt(DHx * DHx + DHy * DHy + DHz * DHz);
                            float HAmag = sqrt(HAx * HAx + HAy * HAy + HAz * HAz);
                            float angle = acos(dot / (DHmag * HAmag)) * 180.0 / M_PI;
                            if (angle > dha_angle) {
                                // get the angle from the donor - acceptor - antecedent using the dot product
                                float DAnx = antecedent->x() - acceptor->x();
                                float DAny = antecedent->y() - acceptor->y();
                                float DAnz = antecedent->z() - acceptor->z();
                                float dot2 = DHx * DAnx + DHy * DAny + DHz * DAnz;
                                float DAnmag = sqrt(DAnx * DAnx + DAny * DAny + DAnz * DAnz);
                                float angle2 = acos(dot2 / (DHmag * DAnmag)) * 180.0 / M_PI;
                                if (angle2 > daa_angle and hydrogen->distance(*acceptor) < distance) {
                                    if (verbose) {
                                        cout << "Hydrogen bond found with:\nDHA: " << angle << " degrees\nDAA: " << angle2 << " degrees\nDistance: " << hydrogen->distance(*acceptor) << " angstroms\nAtoms:\n";
                                        cout << donor->str();
                                        cout << hydrogen->str();
                                        cout << acceptor->str();
                                        cout << antecedent->str() << "\n";
                                    }
                                    PROT::HydrogenBond hbond;
                                    hbond.donor_residue = this;
                                    // of the donor atom is a backbone atom
                                    hbond.donor_backbone = donor->is_backbone_atom();
                                    hbond.acceptor_residue = other;
                                    hbond.acceptor_backbone = acceptor->is_backbone_atom();
                                    hbond.donor_atom = donor;
                                    hbond.hydrogen = hydrogen;
                                    hbond.acceptor_atom = acceptor;
                                    hbond.antecedent = antecedent;
                                    hbond.DHA_angle = angle;
                                    hbond.DAAn_angle = angle2;
                                    hbond.distance = hydrogen->distance(*acceptor);
                                    hbonds.push_back(hbond);
                                    break;
                                }
                            } 
                        } 
                        two_antecedents = false;
                    } 
                }   
            }
        }
    }

    // check the other direction
    for (size_t i_atom = 0; i_atom < other->m_count; i_atom++){
        PROT::Atom* hydrogen = &other->m_atoms[i_atom];
        if (hydrogen->m_element != "H") {
            continue;
        }
        // get the donor strings 
        vector<string> donors = hbond_donors[other_residue];
        // go through donors and if the atom is a donatable hydrogen, get the donor atom
        for (string donor_pair : donors) {
            // split the donor pair, the donor is the second element
            vector<string> hydrogen_and_donor= Text::split(donor_pair, ' ');
            if (Text::contains(hydrogen->name(), hydrogen_and_donor[0])) {
                string donor_name = hydrogen_and_donor[1];
                Text::strip(donor_name);
                // get the donor atom by reference
                PROT::Atom* donor = other->get_atom(donor_name);
                // go through the atoms of the other residue and if they are acceptors, get their acceptor atoms
                for (size_t j_atom = 0; j_atom < m_count; j_atom++){
                    PROT::Atom* acceptor = &m_atoms[j_atom];
                    if (acceptor->m_element != "O" && acceptor->m_element != "N") {
                        continue;
                    }
                    bool two_antecedents = false;
                    // get the acceptor strings
                    vector<string> acceptors = hbond_acceptors[this_residue];
                    // go through acceptors and if the atom is an acceptor, get the acceptor atom
                    for (string acceptor_pair : acceptors) {
                        // split the acceptor pair, the acceptor is the first element
                        vector<string> acceptor_and_antecedent = Text::split(acceptor_pair, ' ');
                        // if the second element has an underscore, then there are 2 antecedents, get their midpoint as the antecedent
                        if (Text::contains(acceptor_and_antecedent[1], "_")) {
                            two_antecedents = true;
                        }
                        if (Text::contains(acceptor->name(), acceptor_and_antecedent[0])) {
                            PROT::Atom dummy;
                            PROT::Atom* antecedent;
                            // if there are two antecedents, get the midpoint
                            if (two_antecedents){
                                string antecedent_name = acceptor_and_antecedent[1];
                                // split the antecedent pair
                                vector<string> antecedents = Text::split(antecedent_name, '_');
                                // get the antecedent atoms by reference
                                Text::strip(antecedents[0]);
                                Text::strip(antecedents[1]);
                                PROT::Atom * antecedent1 = this->get_atom(antecedents[0]);
                                PROT::Atom * antecedent2 = this->get_atom(antecedents[1]);
                                // get the midpoint of the two antecedents as and set them as the antecedent coords
                                dummy.m_coors[0] = (antecedent1->m_coors[0] + antecedent2->m_coors[0]) / 2;
                                dummy.m_coors[1] = (antecedent1->m_coors[1] + antecedent2->m_coors[1]) / 2;
                                dummy.m_coors[2] = (antecedent1->m_coors[2] + antecedent2->m_coors[2]) / 2;
                                dummy.m_name = "dummy";
                                antecedent = &dummy;
                            } else {
                                string antecedent_name = acceptor_and_antecedent[1];
                                Text::strip(antecedent_name);
                                antecedent = this->get_atom(antecedent_name);
                            }
                            // get the angle from the donor - hydrogen - acceptor using the dot product
                            float DHx = donor->x() - hydrogen->x();
                            float DHy = donor->y() - hydrogen->y();
                            float DHz = donor->z() - hydrogen->z();
                            float HAx = acceptor->x() - hydrogen->x();
                            float HAy = acceptor->y() - hydrogen->y();
                            float HAz = acceptor->z() - hydrogen->z();
                            float dot = DHx * HAx + DHy * HAy + DHz * HAz;
                            float DHmag = sqrt(DHx * DHx + DHy * DHy + DHz * DHz);
                            float HAmag = sqrt(HAx * HAx + HAy * HAy + HAz * HAz);
                            float angle = acos(dot / (DHmag * HAmag)) * 180.0 / M_PI;
                            if (angle > dha_angle) {
                                // get the angle from the donor - acceptor - antecedent using the dot product
                                float DAnx = antecedent->x() - acceptor->x();
                                float DAny = antecedent->y() - acceptor->y();
                                float DAnz = antecedent->z() - acceptor->z();
                                float dot2 = DHx * DAnx + DHy * DAny + DHz * DAnz;
                                float DAnmag = sqrt(DAnx * DAnx + DAny * DAny + DAnz * DAnz);
                                float angle2 = acos(dot2 / (DHmag * DAnmag)) * 180.0 / M_PI;
                                if (angle2 > daa_angle and hydrogen->distance(*acceptor) < distance) {
                                    if (verbose) {
                                        cout << "Hydrogen bond found with:\nDHA: " << angle << " degrees\nDAA: " << angle2 << " degrees\nDistance: " << hydrogen->distance(*acceptor) << " angstroms\nAtoms:\n";
                                        cout << donor->str();
                                        cout << hydrogen->str();
                                        cout << acceptor->str();
                                        cout << antecedent->str() << "\n";
                                    }
                                    PROT::HydrogenBond hbond;
                                    hbond.donor_residue = other;
                                    hbond.donor_backbone = donor->is_backbone_atom();
                                    hbond.acceptor_residue = this;
                                    hbond.acceptor_backbone = acceptor->is_backbone_atom();
                                    hbond.donor_atom = donor;
                                    hbond.hydrogen = hydrogen;
                                    hbond.acceptor_atom = acceptor;
                                    hbond.antecedent = antecedent;
                                    hbond.DHA_angle = angle;
                                    hbond.DAAn_angle = angle2;
                                    hbond.distance = hydrogen->distance(*acceptor);
                                    hbonds.push_back(hbond);
                                    break;
                                }
                            }
                        }
                        two_antecedents = false;
                    }
                }
            }
        }
    }

    return hbonds;
}