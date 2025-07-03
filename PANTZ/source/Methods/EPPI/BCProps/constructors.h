/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file implements:
 * A. the constructor of the BCProps class from a vector of pointers to the 
 * interactions and a set of ria values

 * B. constructor of the BCProps class from 
 * a directory containing:
 * 1. features.txt (produced in EPPI/calculate_eppi_features.h)
 * 2. rosetta_residue_scores.sc (produced in EPPI/calculate_eppi_features.h)
 * 3. rosetta_interface_score.sc (produced in EPPI/calculate_eppi_features.h)
  */

// Make sure the file is being included in a compiled program in the expected
// manner
#ifndef EPPI_BCProps_Loading_Status
#error constructors.h must be included by BCProps.h
#endif


// the Constructor from interactions and ria values
void EPPI::BCProps::interaction_constructor(const string& name, vector<EPPI::HydrogenBond*> hydrogen_bonds,
    vector<EPPI::SaltBridge*> salt_bridges, vector<EPPI::Hydrophobic*> hydrophobic_interactions,
    vector<double> ria_values) {
    // set the name
    m_name = name;
    // set the base features
    double n_epp_hb = 0.0; double n_epp_hp = 0.0; double n_epp_sb = 0.0;
    double n_np_hb = 0.0; double n_np_hp = 0.0; double n_np_sb = 0.0;
    double rl_epp_hb = 0.0; double rl_epp_hp = 0.0; double rl_epp_sb = 0.0;
    double rl_np_hb = 0.0; double rl_np_hp = 0.0; double rl_np_sb = 0.0;
    double dg_separated = 0.0; double bsa = 0.0;  double sc = 0.0;
    double epbe_hb = 0.0; double epbe_hp = 0.0; double epbe_sb = 0.0;
    double n_st_st_hb = 0.0; double n_st_st_hp = 0.0; double n_st_st_sb = 0.0;
    // get the ria values
    dg_separated = ria_values[0];
    bsa = ria_values[1];
    sc = ria_values[2];

    for (int i = 0; i < hydrogen_bonds.size(); i++) {
        if (hydrogen_bonds[i]->is_real()) {
            n_epp_hb++;
            rl_epp_hb += (hydrogen_bonds[i]->residue1->pre_free_rotamers - hydrogen_bonds[i]->residue1->bound_free_rotamers) +
                (hydrogen_bonds[i]->residue2->pre_free_rotamers - hydrogen_bonds[i]->residue2->bound_free_rotamers);
            epbe_hb += (0.74*hydrogen_bonds[i]->energy);
            if (hydrogen_bonds[i]->residue1->boundstable && hydrogen_bonds[i]->residue2->boundstable) {
                n_st_st_hb++;
            }
        } else {
            n_np_hb++;
            rl_np_hb += (hydrogen_bonds[i]->residue1->pre_free_rotamers - hydrogen_bonds[i]->residue1->bound_free_rotamers) +
                (hydrogen_bonds[i]->residue2->pre_free_rotamers - hydrogen_bonds[i]->residue2->bound_free_rotamers);
        }
    }
    for (int i = 0; i < salt_bridges.size(); i++) {
        if (salt_bridges[i]->is_real()) {
            n_epp_sb++;
            rl_epp_sb += (salt_bridges[i]->residue1->pre_free_rotamers - salt_bridges[i]->residue1->bound_free_rotamers) +
                (salt_bridges[i]->residue2->pre_free_rotamers - salt_bridges[i]->residue2->bound_free_rotamers);
            epbe_sb += (0.75*salt_bridges[i]->energy);
            if (salt_bridges[i]->residue1->boundstable && salt_bridges[i]->residue2->boundstable) {
                n_st_st_sb++;
            }
        } else {
            n_np_sb++;
            rl_np_sb += (salt_bridges[i]->residue1->pre_free_rotamers - salt_bridges[i]->residue1->bound_free_rotamers) + 
                (salt_bridges[i]->residue2->pre_free_rotamers - salt_bridges[i]->residue2->bound_free_rotamers);
        }
    }
    for (int i = 0; i < hydrophobic_interactions.size(); i++) {
        if (hydrophobic_interactions[i]->is_real()) {
            n_epp_hp++;
            rl_epp_hp += (hydrophobic_interactions[i]->residue1->pre_free_rotamers - hydrophobic_interactions[i]->residue1->bound_free_rotamers) +
                (hydrophobic_interactions[i]->residue2->pre_free_rotamers - hydrophobic_interactions[i]->residue2->bound_free_rotamers);
            epbe_hp += (0.87*hydrophobic_interactions[i]->energy);
            if (hydrophobic_interactions[i]->residue1->boundstable && hydrophobic_interactions[i]->residue2->boundstable) {
                n_st_st_hp++;
            }
        } else {
            n_np_hp++;
            rl_np_hp += (hydrophobic_interactions[i]->residue1->pre_free_rotamers - hydrophobic_interactions[i]->residue1->bound_free_rotamers) +
                (hydrophobic_interactions[i]->residue2->pre_free_rotamers - hydrophobic_interactions[i]->residue2->bound_free_rotamers);
        }
    }

    // set the values
    m_baseFeatures = {n_epp_sb, n_epp_hb, n_epp_hp, 
                       n_np_sb, n_np_hb, n_np_hp, 
                       rl_epp_sb, rl_epp_hb, rl_epp_hp, 
                       rl_np_sb, rl_np_hb, rl_np_hp, 
                       dg_separated, bsa, sc, 
                       epbe_sb, epbe_hb, epbe_hp, 
                       n_st_st_sb, n_st_st_hb, n_st_st_hp};

}

// Implement the constructor function from a name, directory from calculate_eppi_features and a boolean for verbosity
void EPPI::BCProps::directory_constructor (const string& name, const string& dir, const bool verbose) {

    //check that the directory has the 3 necessary files
    vector<string> files = METHODS::listdir(dir);
    if (find(files.begin(), files.end(), "features.txt") == files.end() ||
        find(files.begin(), files.end(), "rosetta_residue_scores.sc") == files.end() ||
        find(files.begin(), files.end(), "rosetta_interface_score.sc") == files.end()) {
        throw PANTZ_error("Error: directory "+dir+" does not have the necessary files");
        return;
    }
    // Parse the per res file
    map<string, double> energy_map = parse_per_res_file(dir+"/rosetta_residue_scores.sc");
    // parse the ria file
    vector<double> ria_values = parse_ria_file(dir+"/rosetta_interface_score.sc");

    // declare the vectors to store the interactions
    vector<EPPI::HydrogenBond*> hydrogen_bonds;
    vector<EPPI::SaltBridge*> salt_bridges;
    vector<EPPI::Hydrophobic*> hydrophobic_interactions;

    // parse the EPPI file
    parse_eppi_file(dir+"/features.txt", hydrogen_bonds, salt_bridges, hydrophobic_interactions);

    // set the energy of the interactions
    for (int i = 0; i < hydrogen_bonds.size(); i++) {
        hydrogen_bonds[i]->set_energy(energy_map);
    }
    for (int i = 0; i < salt_bridges.size(); i++) {
        salt_bridges[i]->set_energy(energy_map);
    }
    for (int i = 0; i < hydrophobic_interactions.size(); i++) {
        hydrophobic_interactions[i]->set_energy(energy_map);
    }
    // create the base features using the interaction constructor
    interaction_constructor(name, hydrogen_bonds, salt_bridges, hydrophobic_interactions, ria_values);
    // If verbose, remove the 3 files
    if (!verbose) {
        string cmd = "rm " + dir + "/features.txt " + dir + "/rosetta_residue_scores.sc " +
                    dir + "/rosetta_interface_score.sc " + dir + "/log.txt " + dir + "/rosetta_output.out";
        system (cmd.c_str());
    }
}