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

using namespace std;
using namespace PROT;

#include <map>
#include <vector>
#include <string>
#include <set>

// // map of the number of rotamers for each amino acid
// map<char, size_t> rotamer_count = {
//     {'A', 1},
//     {'R', 53},
//     {'N', 17},
//     {'D', 9},
//     {'C', 3},
//     {'Q', 31},
//     {'E', 26},
//     {'G', 1},
//     {'H', 9},
//     {'I', 8},
//     {'L', 8},
//     {'K', 47},
//     {'M', 22},
//     {'F', 5},
//     {'P', 0},
//     {'S', 3},
//     {'T', 3},
//     {'W', 9},
//     {'Y', 6},
//     {'V', 3}
// };

// // a function to make rotamers of this residue
// // temporarily just load them from the rotlib
// void PROT::Residue::set_rotamers() {
//     string path = "/home/shared/rjp0029_lab/rotlib/independ/";
//     vector<PROT::Residue> rotamers;
//     // get this 3 letter name in lowercase
//     string name = m_name;
//     Text::lower(name);
//     // if his, change to hsd
//     if (name == "his") {
//         name = "hsd";
//     }
//     // get the 1 letter name
//     char AA1 = this->AA1();
//     // get the number of rotamers for this amino acid
//     size_t count = rotamer_count[AA1];
//     // go through each and align, and add to list
//     for (size_t i = 0; i < count; i++) {
//         // get the rotamer name
//         string rotamer_filename = path + name + to_string(i+1) + ".pdb";
//         // ifstream to read the rotamer file
//         ifstream rotamer_file(rotamer_filename);
//         vector<PROT::Atom*> rotamer_atoms;
//         for (string line; getline(rotamer_file, line);) {
//             // if not an atom line, skip it
//             if (line.substr(0, 4) != "ATOM") {
//                 continue;
//             }
//             // create a new atom from the line
//             PROT::Atom* atom = new PROT::Atom(line);
//             // set the m_element to the first character of the atom name
//             atom->m_element = atom->m_name[0];
//             // add it to the residue
//             rotamer_atoms.push_back(atom);
//         } 
//         // create a residue from the atoms
//         PROT::Residue rotamer(rotamer_atoms);
//         // align the rotamer to this residue by moving to origin and aligning N, CB
//         PROT::Matrix this_residue_mover(1, 3);
//         this->position(this_residue_mover, true);
//         PROT::Matrix rotamer_mover(1, 3);
//         rotamer.position(rotamer_mover, true);
//         // move everything back to the original position
//         this->position(this_residue_mover, false);
//         rotamer.position(this_residue_mover, false);
//         // add it to the list
//         rotamers.push_back(rotamer);
//     }
//     m_rotamers = rotamers;
// }

map<char, vector<vector<string>>> chi_definitions = {
    {'A', {}},
    {'R', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD"}, {"CB", "CG", "CD", "NE"}, {"CG", "CD", "NE", "CZ"}}},
    {'N', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "OD1"}}},
    {'D', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "OD1"}}},
    {'Q', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD"}, {"CB", "CG", "CD", "OE1"}}},
    {'E', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD"}, {"CB", "CG", "CD", "OE1"}}},
    {'G', {}},
    {'H', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "ND1"}}},
    {'K', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD"}, {"CB", "CG", "CD", "CE"}, {"CG", "CD", "CE", "NZ"}}},
    {'M', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "SD"}, {"CB", "CG", "SD", "CE"}}},
    {'F', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD1"}}},
    {'W', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD1"}}},
    {'Y', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD1"}}},
    {'I', {{"N", "CA", "CB", "CG1"}, {"CA", "CB", "CG1", "CD"}}},
    {'L', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD1"}}},
    {'P', {{"N", "CA", "CB", "CG"}, {"CA", "CB", "CG", "CD"}, {"CB", "CG", "CD", "HD2"}}},
    {'S', {{"N", "CA", "CB", "OG"}}},
    {'T', {{"N", "CA", "CB", "OG1"}}},
    {'V', {{"N", "CA", "CB", "CG1"}}},
    {'C', {{"N", "CA", "CB", "SG"}}}
};

// average phi and psi for 10000 residues from 1000 random pdb files
map<char, vector<float>> default_dihedrals = {
    {'A', {-80, 30}},
    {'R', {-80, 30}},
    {'N', {-70, 30}},
    {'D', {-80, 30}},
    {'Q', {-80, 30}},
    {'E', {-80, 20}},
    {'G', {-80, 30}},
    {'H', {-90, 40}},
    {'I', {-90, 50}},
    {'L', {-80, 30}},
    {'K', {-80, 30}},
    {'M', {-80, 30}},
    {'F', {-90, 50}},
    {'P', {-70, 70}},
    {'S', {-90, 50}},
    {'T', {-90, 60}},
    {'W', {-90, 40}},
    {'Y', {-90, 50}},
    {'V', {-90, 60}},
    {'C', {-90, 60}}
};

// function to get current chi angles
vector<float> get_chi_angles(PROT::Residue * res, bool verbose) {
    vector<float> chi_angles;
    // get the chi angles
    char AA1 = res->AA1();
    // chi angle atom definitions
    PROT::Atom* atom1;
    PROT::Atom* atom2;
    PROT::Atom* atom3;
    PROT::Atom* atom4;
    for (size_t i = 0; i < chi_definitions[AA1].size(); i++) {
        // get the atom names
        vector<string> atom_names = chi_definitions[AA1][i];
        // get the atoms
        atom1 = res->get_atom(atom_names[0]);
        atom2 = res->get_atom(atom_names[1]);
        atom3 = res->get_atom(atom_names[2]);
        atom4 = res->get_atom(atom_names[3]);
        // calculate the dihedral angle
        float chi = calculate_dihedral(*atom1, *atom2, *atom3, *atom4);
        // append the chi angle to the list
        chi_angles.push_back(chi);
        if (verbose) {
            cout<<"\tChi angle "<<i+1<<": "<<chi<<"\t";
        }
    }
    return chi_angles;
}

// get the atom index in a residue
size_t get_atom_index(PROT::Residue * res, string name) {
    for (size_t i = 0; i < res->size(); i++) {
        if (res->get_atom(i)->name() == name) {
            return i;
        }
    }
    return -999;
}

// the reference level heirarchy of the atoms (A, B, G, D, E, etc.)
string level_heirarchy = "ABGDEZ";

// a function to get the atoms before a given chi angle
set<size_t> get_fixed_atoms(PROT::Residue * res, size_t chi_index) {
    set<size_t> fixed_atoms;
    // get the chi angles
    char AA1 = res->AA1();
    // chi angle atom definitions
    PROT::Atom* atom1;
    PROT::Atom* atom2;
    PROT::Atom* atom3;
    PROT::Atom* atom4;
    for (size_t i = 0; i < chi_index+1; i++) {
        // get the atom names
        vector<string> atom_names = chi_definitions[AA1][i];
        // get the atoms
        atom1 = res->get_atom(atom_names[0]);
        atom2 = res->get_atom(atom_names[1]);
        // get the level of the second atom (Beta, Gamma, Delta, Epsilon, etc.)
        // since none of the second atoms have numbers in their names, we can just get the second character
        string level = atom_names[1].substr(1, 1);
        // cout<<i<<" level: "<<level<<endl;
        // go through the atoms in the residue and fix all atoms up to the level of the second atom (thoose with the same second character or less)
        for (size_t j = 0; j < res->size(); j++) {
            // get the atom name
            string atom_name = res->get_atom(j)->name();
            // get the level of the atom
            string atom_level = atom_name.substr(1, 1);
            // cout<<"atom name: "<<atom_name<<" atom level: "<<atom_level<<endl;
            // if the atom level is less than or equal to the second atom level, add it to the fixed atoms
            if (level_heirarchy.find(atom_level) <= level_heirarchy.find(level)) {
                fixed_atoms.insert(j);
            }
        }
    }
    return fixed_atoms;
}

// a function to load a residue from a file with a single residue in it
PROT::Residue load_residue(string filename) {
    ifstream file(filename);
    vector<PROT::Atom*> atoms;
    for (string line; getline(file, line);) {
        // if not an atom line, skip it
        if (line.substr(0, 4) != "ATOM") {
            continue;
        }
        // create a new atom from the line
        PROT::Atom* atom = new PROT::Atom(line);
        // add it to the residue
        atoms.push_back(atom);
    }
    // create a residue from the atoms
    PROT::Residue residue(atoms);
    return residue;
}

void PROT::Residue::set_rotamers() {
    // ensure the phi and psi angles of this residue are set
    if (m_phi < -999 || m_psi < -999) {
        string error = "Phi and psi angles must be set before rotamers can be "
                       "generated, call calculate_dihedrals() on the Protein\n";
        // use the default dihedrals for this residue
        m_phi = default_dihedrals[AA1()][0];
        m_psi = default_dihedrals[AA1()][1];
    }
    // get the 3 letter name of this residue
    string name = m_name;
    // lowercase the name
    Text::lower(name);
    // if the residue is alanine or glycine, return
    if (name == "ala") {
        // const string path = string(ROTLIB_PATH) + "/" + name + "_avg_phi_psi_rotamer.pdb";
        // // load the default rotamers and then return
        // PROT::Residue residue = load_residue(path);
        // // change the number of the residue to this residues number
        // residue.set_number(m_number, ' ', false);
        // residue.set_number(m_internal, ' ', true);
        // residue.set_protein(m_protein);
        vector<PROT::Residue> rotamers;
        // rotamers.push_back(residue);
        // add this residue
        rotamers.push_back(*this);
        m_rotamers = rotamers;
        return;
    } else if (name == "gly") {
        // const string path = string(ROTLIB_PATH) + "/" + name + "_avg_phi_psi_rotamer.pdb";
        // cout<<path<<endl;
        // // load the default rotamers and then return
        // PROT::Residue residue = load_residue(path);
        // cout<<"loaded residue"<<endl;
        // // change the number of the residue to this residues number
        // residue.set_number(m_number, ' ', false);
        // residue.set_number(m_internal, ' ', true);
        // residue.set_protein(m_protein);
        // cout<<residue.str()<<endl;
        vector<PROT::Residue> rotamers;
        // rotamers.push_back(residue);
        // add this residue
        rotamers.push_back(*this);
        m_rotamers = rotamers;
        return;
    } else if (name == "pro") {
        // const string path = string(ROTLIB_PATH) + "/" + name + "_avg_phi_psi_rotamer.pdb";
        // // load the default rotamers and then return
        // PROT::Residue residue = load_residue(path);
        // // change the number of the residue to this residues number
        // residue.set_number(m_number, ' ', false);
        // residue.set_number(m_internal, ' ', true);
        // residue.set_protein(m_protein);
        vector<PROT::Residue> rotamers;
        // rotamers.push_back(residue);
        rotamers.push_back(*this);
        m_rotamers = rotamers;
        return;
    }
    // otherwise load the rotamer library file
    ifstream rotlib(string(PANTZ_PATH) + "/source/external/rotamer_library/ExtendedOpt1-5/" + name + ".bbdep.rotamers.lib");
    // if the file is not open, throw an error
    if (!rotlib.is_open()) {
        string error = "Could not open rotamer library file for residue " + name + "\n";
        throw PANTZ_error (
            error
        );
    }

    // get phi and psi agles rounded to nearest ten
    int phi = round(m_phi/10)*10;
    int psi = round(m_psi/10)*10;

    // chi angle atom definitions
    PROT::Atom* atom1;
    PROT::Atom* atom2;
    PROT::Atom* atom3;
    PROT::Atom* atom4;
    // get the chi angles
    vector<float> current_chi_angles;
    char AA1 = this->AA1();

    for (size_t i = 0; i < chi_definitions[AA1].size(); i++) {
        // get the atom names
        vector<string> atom_names = chi_definitions[AA1][i];
        // get the atoms
        atom1 = get_atom(atom_names[0]);
        atom2 = get_atom(atom_names[1]);
        atom3 = get_atom(atom_names[2]);
        atom4 = get_atom(atom_names[3]);
        // calculate the dihedral angle
        float chi = calculate_dihedral(*atom1, *atom2, *atom3, *atom4);
        // append the chi angle to the list
        current_chi_angles.push_back(chi);
    }

    vector<vector<float>> all_chi_angles;
    // print the lines of this file
    for (string line; getline(rotlib, line);) {
        if (line[0] == '#') {
            continue;
        }
        string phi_s = line.substr(5, 5);
        string psi_s = line.substr(10, 5);
        // if the second and third element of the line are the phi and psi angles print the line
        if (stoi(phi_s) == phi && stoi(psi_s) == psi) {
            // get the chi angles
            vector<float> chis {stof(line.substr(47, 6)), stof(line.substr(55, 6)), stof(line.substr(63, 6)), stof(line.substr(71, 6))};
            all_chi_angles.push_back(chis);
        }
    }

    // reduce trailing 0's from the chi angles
    vector<vector<float>> all_chi_angles_reduced;
    for (vector<float> chi_angles : all_chi_angles) {
        vector<float> chi_angles_reduced;
        for (float chi : chi_angles) {
            if (chi != 0) {
                chi_angles_reduced.push_back(chi);
            }
        }
        all_chi_angles_reduced.push_back(chi_angles_reduced);
    }

    // Create a vector to store Rotamer objects
    vector<PROT::Residue> rotamers;
    vector<float> new_chi_angles;
    // Iterate through each set of chi angles
    for (size_t r = 0; r < all_chi_angles_reduced.size(); ++r) {
        // Get the chi angles for this rotamer
        vector<float> chi_angles = all_chi_angles_reduced[r];
        // Create a new Rotamer object that is a copy of the original residue
        PROT::Residue rotamer = *this;
        // get the current chi angles
        current_chi_angles = get_chi_angles(&rotamer, false);
        // indices of the atoms to fix
        set<size_t> fixed_atoms;
        // add backbone atoms to fixed atoms
        for (size_t i = 0; i < m_count; i++) {
            if (rotamer.m_atoms[i].is_backbone_atom()) {
                fixed_atoms.insert(i);
            }
        }
        // Apply rotations to the sidechain atoms based on chi angles
        for (size_t i = 0; i < chi_angles.size(); ++i) {
            // atoms to rotate around
            vector<string> atom_names = chi_definitions[AA1][i];
            // get the atoms
            atom1 = rotamer.get_atom(atom_names[1]);
            atom2 = rotamer.get_atom(atom_names[2]);
            // Add atoms before the current rotatable bond to the fixed list
            set<size_t> new_fixed_atoms = get_fixed_atoms(&rotamer, i);
            fixed_atoms.insert(new_fixed_atoms.begin(), new_fixed_atoms.end());
            // // angle is the chi angle in radians
            const coor angle = (chi_angles[i] - current_chi_angles[i]) * M_PI / 180.0f;
            // vector is the bond vector of the atom to rotate around
            const coor vector[3] = {atom2->m_coors[0] - atom1->m_coors[0], atom2->m_coors[1] - atom1->m_coors[1], atom2->m_coors[2] - atom1->m_coors[2]};
            // normalize the vector
            coor norm = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
            const coor vector_norm[3] = {vector[0]/norm, vector[1]/norm, vector[2]/norm};
            // // Create a matrix for the rotation
            PROT::Matrix rotation_matrix(angle, vector_norm);
            // Apply rotation to all atoms after the rotation center
            for (size_t j = 0; j < m_count; ++j) {
                // dont rotate fixed atoms
                if (find(fixed_atoms.begin(), fixed_atoms.end(), j) != fixed_atoms.end()) {
                    // cout<<"fixed atom for chi angle "<<i+1<<": "<<rotamer.m_atoms[j].name()<<endl;
                    continue;
                }
                // translate the atom coordinates to the origin
                rotamer.m_atoms[j].m_coors[0] -= atom1->m_coors[0];
                rotamer.m_atoms[j].m_coors[1] -= atom1->m_coors[1];
                rotamer.m_atoms[j].m_coors[2] -= atom1->m_coors[2];
                // apply the rotation
                rotamer.m_atoms[j].rotate(rotation_matrix);
                // translate the atom coordinates back to the original position
                rotamer.m_atoms[j].m_coors[0] += atom1->m_coors[0];
                rotamer.m_atoms[j].m_coors[1] += atom1->m_coors[1];
                rotamer.m_atoms[j].m_coors[2] += atom1->m_coors[2];
            }
            // reset current chi angles
            current_chi_angles = get_chi_angles(&rotamer, false);
        }
        // Add the generated rotamer to the lis t
        rotamers.push_back(rotamer);
        // cout<<"final chi angles for rotamer "<<r<<": ";
        // get_chi_angles(&rotamer, true);
        // cout<<endl;
        // break;
    }
    // Set the rotamers of this residue to the generated rotamers
    m_rotamers = rotamers;
}
