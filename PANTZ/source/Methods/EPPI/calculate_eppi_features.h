/* Created by clay at Auburn University.
*
* This file implements the function to calculate the expected persistent pairwise
* interaction features between two proteins. This is based on the paper by
* Richard et al. 2024 (https://doi.org/10.1002/prot.26773)
*
*/

// This file is supposed to be loaded by EPPI.h
#ifndef EPPI_Loading_Status
#error EPPI methods must be included by EPPI.h
#endif

// a function to scan the interface and get the expected persistent pairwise interaction features
void EPPI::calculate_eppi_features(vector<PROT::Protein*>& proteins, string& interface, string& output_path, bool verbose){
    // update atoms after rosetta
    for (size_t i = 0; i < proteins.size(); i++){
        proteins[i]->update_atoms_after_Rosetta();
    }
    // get the characters before the _ in the interface string
    string interface_side1 = interface.substr(0, interface.find("_"));
    string interface_side2 = interface.substr(interface.find("_")+1, interface.size());
    // get the interface atoms for the kdtree
    vector<PROT::Atom*> interface_atoms_ptrs;
    vector<PROT::Residue*> protein1_residues;
    vector<PROT::Residue*> protein2_residues;
    float cutoff = 6.0; // cutoff for interface residues
    // get the residues in interface side 1 that are within cutoff of a residue in interface side 2
    for (size_t i = 0; i < proteins.size(); i++){
        // go through residues in the protein
        for (size_t j = 0; j < proteins[i]->size(); j++){
            PROT::Residue * res = proteins[i]->operator()(j, ' ', true);
            // ensure they are on side 1
            if (interface_side1.find(res->protein()) != string::npos){
                // go through residues in the protein again
                for (size_t k = 0; k < proteins.size(); k++){
                    // ensure this protein is on side 2
                    if (interface_side2.find(proteins[k]->name()) != string::npos){
                        // go through the residues in the protein
                        for (size_t l = 0; l < proteins[k]->size(); l++){
                            PROT::Residue * res2 = proteins[k]->operator()(l, ' ', true);
                            if (res->min_distance(*res2) < cutoff){
                                protein1_residues.push_back(proteins[i]->operator()(j, ' ', true));
                                for (size_t m = 0; m < res->size(); m++){
                                    interface_atoms_ptrs.push_back(res->get_atom(m));
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // get the residues in interface side 2 that are within cutoff of a residue in interface side 1
    for (size_t i = 0; i < proteins.size(); i++){
        // go through residues in the protein
        for (size_t j = 0; j < proteins[i]->size(); j++){
            PROT::Residue * res = proteins[i]->operator()(j, ' ', true);
            // ensure they are on side 2
            if (interface_side2.find(res->protein()) != string::npos){
                // go through residues in the protein again
                for (size_t k = 0; k < proteins.size(); k++){
                    // ensure this protein is on side 1
                    if (interface_side1.find(proteins[k]->name()) != string::npos){
                        // go through the residues in the protein
                        for (size_t l = 0; l < proteins[k]->size(); l++){
                            PROT::Residue * res2 = proteins[k]->operator()(l, ' ', true);
                            if (res->min_distance(*res2) < cutoff){
                                protein2_residues.push_back(res);
                                for (size_t m = 0; m < res->size(); m++){
                                    PROT::Atom * atom = res->get_atom(m);
                                    interface_atoms_ptrs.push_back(atom);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // remove duplicates in the residue lists
    sort(protein1_residues.begin(), protein1_residues.end());
    protein1_residues.erase(unique(protein1_residues.begin(), protein1_residues.end()), protein1_residues.end());
    sort(protein2_residues.begin(), protein2_residues.end());
    protein2_residues.erase(unique(protein2_residues.begin(), protein2_residues.end()), protein2_residues.end());

    // remove duplicates in the atoms
    sort(interface_atoms_ptrs.begin(), interface_atoms_ptrs.end());
    interface_atoms_ptrs.erase(unique(interface_atoms_ptrs.begin(), interface_atoms_ptrs.end()), interface_atoms_ptrs.end());

    // set the sasa points for the interface atoms
    METHODS::set_sasa_points(interface_atoms_ptrs);
    
    // intialize the interactions
    vector<PROT::HydrogenBond> hbonds;
    vector<PROT::SaltBridge> salt_bridges;
    vector<PROT::Hydrophobic> hydrophobic_interactions;
    // the values used in clays paper
    float sb_distance = 4.0; // distance threshold for salt bridge in angstroms
    float hb_distance = 2.5; // distance threshold for hydrogen bond in angstroms
    float dha_angle = 120; // angle threshold for hydrogen bond in degrees
    float daa_angle = 90; // angle threshold for hydrogen bond in degrees
    float sasa_cutoff = 24.438; // cutoff for hydrophobic interactions
    for (auto res1 : protein1_residues){
        for (auto res2 : protein2_residues){
            // get the number of hydrogen bonds between the residues
            vector<PROT::HydrogenBond> hbond = res1->hbond(res2, hb_distance, dha_angle, daa_angle, verbose);
            vector<PROT::SaltBridge> salt_bridge = res1->salt_bridge(res2, sb_distance, verbose);
            vector<PROT::Hydrophobic> hydrophobic_interaction = res1->hydrophobic(res2, sasa_cutoff, verbose);  
            // add the hydrogen bonds to the list
            hbonds.insert(hbonds.end(), hbond.begin(), hbond.end());
            salt_bridges.insert(salt_bridges.end(), salt_bridge.begin(), salt_bridge.end());
            hydrophobic_interactions.insert(hydrophobic_interactions.end(), hydrophobic_interaction.begin(), hydrophobic_interaction.end());
        }
    }

    if (verbose){
        cout<<"Number of hydrogen bonds: "<<hbonds.size()<<endl;
        cout<<"Number of salt bridges: "<<salt_bridges.size()<<endl;
        cout<<"Number of hydrophobic interactions: "<<hydrophobic_interactions.size()<<endl;
        cout<<"\nSetting the stability of the residues"<<endl;
    }
    string how;
    bool backbone = false;
    vector<PROT::Residue*> all_residues;
    all_residues.insert(all_residues.end(), protein1_residues.begin(), protein1_residues.end());
    all_residues.insert(all_residues.end(), protein2_residues.begin(), protein2_residues.end());
    // go through the hydrogen bond, salt bridge, and hydrophobic interaction lists and set the stability of the residues
    for (size_t i = 0; i < hbonds.size(); i++){
        // get the neighbors of the donor and acceptor, inter and intra
        // if the donor name is the same as prot1 pass prot1 as the protein, else pass prot2
        vector<PROT::Residue*> donor_intra_neighbors;
        vector<PROT::Residue*> acceptor_intra_neighbors;
        vector<PROT::Residue*> donor_inter_neighbors;
        vector<PROT::Residue*> acceptor_inter_neighbors;
        // if (hbonds[i].donor_residue->protein() == prot1->name()){
        // if tthe donor name char is found in the first chain interface side
        if (interface_side1.find(hbonds[i].donor_residue->protein()) != string::npos){
            donor_intra_neighbors = hbonds[i].donor_residue->get_intra_neighbors(protein1_residues, 14);
            acceptor_intra_neighbors = hbonds[i].acceptor_residue->get_intra_neighbors(protein2_residues, 14);
            donor_inter_neighbors = hbonds[i].donor_residue->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = hbonds[i].acceptor_residue->get_inter_neighbors_res(all_residues, 14);
        } else {
            donor_intra_neighbors = hbonds[i].donor_residue->get_intra_neighbors(protein2_residues, 14);
            acceptor_intra_neighbors = hbonds[i].acceptor_residue->get_intra_neighbors(protein1_residues, 14);
            donor_inter_neighbors = hbonds[i].donor_residue->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = hbonds[i].acceptor_residue->get_inter_neighbors_res(all_residues, 14);
        }
        // get the acceptor rotamers
        vector<PROT::Residue> acceptor_rotamers = hbonds[i].acceptor_residue->get_rotamers();
        vector<PROT::Residue> donor_rotamers = hbonds[i].donor_residue->get_rotamers();
        // set the stability of the residues
        hbonds[i].donor_residue->set_stability(donor_rotamers, donor_intra_neighbors, backbone=hbonds[i].donor_backbone, how="pre");
        hbonds[i].donor_residue->set_stability(donor_rotamers, donor_inter_neighbors, backbone=hbonds[i].donor_backbone, how="bound");
        hbonds[i].acceptor_residue->set_stability(acceptor_rotamers, acceptor_intra_neighbors, backbone=hbonds[i].acceptor_backbone, how="pre");
        hbonds[i].acceptor_residue->set_stability(acceptor_rotamers, acceptor_inter_neighbors, backbone=hbonds[i].acceptor_backbone, how="bound");
        // set the free rotamers
        hbonds[i].donor_free_rot_pre = hbonds[i].donor_residue->free_rotamers(donor_rotamers, donor_intra_neighbors);
        hbonds[i].donor_free_rot_bound = hbonds[i].donor_residue->free_rotamers(donor_rotamers, donor_inter_neighbors);
        hbonds[i].acceptor_free_rot_pre = hbonds[i].acceptor_residue->free_rotamers(acceptor_rotamers, acceptor_intra_neighbors);
        hbonds[i].acceptor_free_rot_bound = hbonds[i].acceptor_residue->free_rotamers(acceptor_rotamers, acceptor_inter_neighbors);
    }
    // set the stability of the salt bridges
    for (size_t i = 0; i < salt_bridges.size(); i++){
        // get the neighbors of the donor and acceptor, inter and intra
        // if the donor name is the same as prot1 pass prot1 as the protein, else pass prot2
        vector<PROT::Residue*> donor_intra_neighbors;
        vector<PROT::Residue*> acceptor_intra_neighbors;
        vector<PROT::Residue*> donor_inter_neighbors;
        vector<PROT::Residue*> acceptor_inter_neighbors;
        // if the donor name char is found in the first chain interface side
        if (interface_side1.find(salt_bridges[i].donor_residue->protein()) != string::npos){
            donor_intra_neighbors = salt_bridges[i].donor_residue->get_intra_neighbors(protein1_residues, 14);
            acceptor_intra_neighbors = salt_bridges[i].acceptor_residue->get_intra_neighbors(protein2_residues, 14);
            donor_inter_neighbors = salt_bridges[i].donor_residue->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = salt_bridges[i].acceptor_residue->get_inter_neighbors_res(all_residues, 14);
        } else {
            donor_intra_neighbors = salt_bridges[i].donor_residue->get_intra_neighbors(protein2_residues, 14);
            acceptor_intra_neighbors = salt_bridges[i].acceptor_residue->get_intra_neighbors(protein1_residues, 14);
            donor_inter_neighbors = salt_bridges[i].donor_residue->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = salt_bridges[i].acceptor_residue->get_inter_neighbors_res(all_residues, 14);
        }
        // get the acceptor rotamers
        vector<PROT::Residue> acceptor_rotamers = salt_bridges[i].acceptor_residue->get_rotamers();
        vector<PROT::Residue> donor_rotamers = salt_bridges[i].donor_residue->get_rotamers();
        // set the stability of the residues
        salt_bridges[i].donor_residue->set_stability(donor_rotamers, donor_intra_neighbors, backbone=salt_bridges[i].donor_backbone, how="pre");
        salt_bridges[i].donor_residue->set_stability(donor_rotamers, donor_inter_neighbors, backbone=salt_bridges[i].donor_backbone, how="bound");
        salt_bridges[i].acceptor_residue->set_stability(acceptor_rotamers, acceptor_intra_neighbors, backbone=salt_bridges[i].acceptor_backbone, how="pre");
        salt_bridges[i].acceptor_residue->set_stability(acceptor_rotamers, acceptor_inter_neighbors, backbone=salt_bridges[i].acceptor_backbone, how="bound");
        // set the free rotamers
        salt_bridges[i].donor_free_rot_pre = salt_bridges[i].donor_residue->free_rotamers(donor_rotamers, donor_intra_neighbors);
        salt_bridges[i].donor_free_rot_bound = salt_bridges[i].donor_residue->free_rotamers(donor_rotamers, donor_inter_neighbors);
        salt_bridges[i].acceptor_free_rot_pre = salt_bridges[i].acceptor_residue->free_rotamers(acceptor_rotamers, acceptor_intra_neighbors);
        salt_bridges[i].acceptor_free_rot_bound = salt_bridges[i].acceptor_residue->free_rotamers(acceptor_rotamers, acceptor_inter_neighbors);
    }
    // set the stability of the hydrophobic interactions
    for (size_t i = 0; i < hydrophobic_interactions.size(); i++){
        // get the neighbors of the donor and acceptor, inter and intra
        // if the donor name is the same as prot1 pass prot1 as the protein, else pass prot2
        vector<PROT::Residue*> donor_intra_neighbors;
        vector<PROT::Residue*> acceptor_intra_neighbors;
        vector<PROT::Residue*> donor_inter_neighbors;
        vector<PROT::Residue*> acceptor_inter_neighbors;
        // if (hbonds[i].donor_residue->protein() == prot1->name()){
        // if tthe donor name char is found in the first chain interface side
        if (interface_side1.find(hydrophobic_interactions[i].res1->protein()) != string::npos){
            donor_intra_neighbors = hydrophobic_interactions[i].res1->get_intra_neighbors(protein1_residues, 14);
            acceptor_intra_neighbors = hydrophobic_interactions[i].res2->get_intra_neighbors(protein2_residues, 14);
            donor_inter_neighbors = hydrophobic_interactions[i].res1->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = hydrophobic_interactions[i].res2->get_inter_neighbors_res(all_residues, 14);
        } else {
            donor_intra_neighbors = hydrophobic_interactions[i].res1->get_intra_neighbors(protein2_residues, 14);
            acceptor_intra_neighbors = hydrophobic_interactions[i].res2->get_intra_neighbors(protein1_residues, 14);
            donor_inter_neighbors = hydrophobic_interactions[i].res1->get_inter_neighbors_res(all_residues, 14);
            acceptor_inter_neighbors = hydrophobic_interactions[i].res2->get_inter_neighbors_res(all_residues, 14);
        }
        // get the acceptor rotamers
        vector<PROT::Residue> acceptor_rotamers = hydrophobic_interactions[i].res2->get_rotamers();
        vector<PROT::Residue> donor_rotamers = hydrophobic_interactions[i].res1->get_rotamers();
        // set the stability of the residues
        hydrophobic_interactions[i].res1->set_stability(donor_rotamers, donor_intra_neighbors, backbone=hydrophobic_interactions[i].res1_backbone, how="pre");
        hydrophobic_interactions[i].res1->set_stability(donor_rotamers, donor_inter_neighbors, backbone=hydrophobic_interactions[i].res1_backbone, how="bound");
        hydrophobic_interactions[i].res2->set_stability(acceptor_rotamers, acceptor_intra_neighbors, backbone=hydrophobic_interactions[i].res2_backbone, how="pre");
        hydrophobic_interactions[i].res2->set_stability(acceptor_rotamers, acceptor_inter_neighbors, backbone=hydrophobic_interactions[i].res2_backbone, how="bound");
        // set the free rotamers
        hydrophobic_interactions[i].res1_free_rot_pre = hydrophobic_interactions[i].res1->free_rotamers(donor_rotamers, donor_intra_neighbors);
        hydrophobic_interactions[i].res1_free_rot_bound = hydrophobic_interactions[i].res1->free_rotamers(donor_rotamers, donor_inter_neighbors);
        hydrophobic_interactions[i].res2_free_rot_pre = hydrophobic_interactions[i].res2->free_rotamers(acceptor_rotamers, acceptor_intra_neighbors);
        hydrophobic_interactions[i].res2_free_rot_bound = hydrophobic_interactions[i].res2->free_rotamers(acceptor_rotamers, acceptor_inter_neighbors);
    }

    // make the output directory if it doesnt exist
    system(("mkdir -p " + output_path).c_str());

    // write the features to a file
    ofstream features_file(output_path + "/features.txt");
    
    for (size_t i = 0; i < hbonds.size(); i++){
        features_file<<"Hydrogen Bond ";
        features_file<<hbonds[i].donor_residue->name()<<"_"<<hbonds[i].donor_residue->number()<<"_"<<hbonds[i].donor_residue->protein()<< " ";
        // if prestable
        if (hbonds[i].donor_residue->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (hbonds[i].donor_residue->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (hbonds[i].donor_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<hbonds[i].donor_free_rot_pre<<" ";
        features_file<<hbonds[i].donor_free_rot_bound<<" ";
        features_file<<", ";
        features_file<<hbonds[i].acceptor_residue->name()<<"_"<<hbonds[i].acceptor_residue->number()<<"_"<<hbonds[i].acceptor_residue->protein()<< " ";
        if (hbonds[i].acceptor_residue->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (hbonds[i].acceptor_residue->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (hbonds[i].acceptor_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<hbonds[i].acceptor_free_rot_pre<<" ";
        features_file<<hbonds[i].acceptor_free_rot_bound<<" ";
        features_file<<", ";
        // output the hbond angle and distance
        features_file<<hbonds[i].distance<<" "<<hbonds[i].DHA_angle<<" "<<hbonds[i].DAAn_angle<<endl;
    }
    
    for (size_t i = 0; i < salt_bridges.size(); i++){
        features_file<<"Salt Bridge ";
        features_file<<salt_bridges[i].donor_residue->name()<<"_"<<salt_bridges[i].donor_residue->number()<<"_"<<salt_bridges[i].donor_residue->protein()<< " ";
        // if prestable
        if (salt_bridges[i].donor_residue->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (salt_bridges[i].donor_residue->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (salt_bridges[i].donor_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<salt_bridges[i].donor_free_rot_pre<<" ";
        features_file<<salt_bridges[i].donor_free_rot_bound<<" ";
        features_file<<", ";
        features_file<<salt_bridges[i].acceptor_residue->name()<<"_"<<salt_bridges[i].acceptor_residue->number()<<"_"<<salt_bridges[i].acceptor_residue->protein()<< " ";
        if (salt_bridges[i].acceptor_residue->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (salt_bridges[i].acceptor_residue->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (salt_bridges[i].acceptor_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<salt_bridges[i].acceptor_free_rot_pre<<" ";
        features_file<<salt_bridges[i].acceptor_free_rot_bound<<" ";
        // output the distance between charged atoms
        features_file<<", "<< salt_bridges[i].distance<<endl;
    }
    
    for (size_t i = 0; i < hydrophobic_interactions.size(); i++){
        features_file<<"Hydrophobic Interaction ";
        features_file<<hydrophobic_interactions[i].res1->name()<<"_"<<hydrophobic_interactions[i].res1->number()<<"_"<<hydrophobic_interactions[i].res1->protein()<< " ";
        // if prestable
        if (hydrophobic_interactions[i].res1->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (hydrophobic_interactions[i].res1->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (hydrophobic_interactions[i].res1_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<hydrophobic_interactions[i].res1_free_rot_pre<<" ";
        features_file<<hydrophobic_interactions[i].res1_free_rot_bound<<" ";
        features_file<<", ";
        features_file<<hydrophobic_interactions[i].res2->name()<<"_"<<hydrophobic_interactions[i].res2->number()<<"_"<<hydrophobic_interactions[i].res2->protein()<< " ";
        if (hydrophobic_interactions[i].res2->prestable()){
            features_file<<"prestable ";
        } else {
            features_file<<"not_prestable ";
        }
        // if bound_stable
        if (hydrophobic_interactions[i].res2->bound_stable()){
            features_file<<"bound_stable ";
        } else {
            features_file<<"not_bound_stable ";
        }
        if (hydrophobic_interactions[i].res2_backbone){
            features_file<<"backbone ";
        } else {
            features_file<<"sidechain ";
        }
        // output the free rotamers
        features_file<<hydrophobic_interactions[i].res2_free_rot_pre<<" ";
        features_file<<hydrophobic_interactions[i].res2_free_rot_bound<<" ";
        // output the bsa
        features_file<<", "<< hydrophobic_interactions[i].phobic_BSASA<<endl;
    }
    features_file.close();
    if (verbose) {
        cout<<"Features written to features.txt"<<endl;
    }

    // run the per residue calculations
    // insert a space between the interface characters break on the _ character for ria
    string interface_;
    for (size_t i = 0; i < interface.size(); i++) {
        // if the character is _ break 
        if (interface[i] == '_') {
            break;
        } else {
            interface_ += interface[i];
            if (i != interface.size() - 1) {
                interface_ += " ";
            }
        }
    }
    // get the proteins 
    vector<PROT::Protein> proteins_;
    for (size_t i = 0; i < proteins.size(); i++) {
        proteins_.push_back(*proteins[i]);
    }
    // run per residue
    if (verbose) {
        cout<<"Running per residue analysis on the original structure"<<endl;
    }
    ofstream output(output_path+"/log.txt");
    Rosetta::Per_Residue(proteins_, output, output_path);
    // run RIA
    if (verbose) {
        cout<<"Running interface analysis on the original structure"<<endl;
    }
    Rosetta::Interface_Analyzer(proteins_, output, output_path, "interface analysis " + interface_);

}