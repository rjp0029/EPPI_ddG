/* Created by clay at Auburn University.
 *
 * This file implements the PRediction of ddG values upon mutation 
 * using the Expected Persistent Pairwise Interaction Features from 
 * Clay's 2025 paper. 
 * 
*/

// This file is supposed to be included by Protocol.h
#ifndef Protocol_Loading_Status
#error General Protocol functions must be included by Protocol.h
#endif

// Function to check if output files already exist
bool already_calculated(string output_path) {
    // if the directory doesnt exist, return false
    ifstream file(output_path);
    if (!file) return false;
    vector<string> files = METHODS::listdir(output_path);
    return (find(files.begin(), files.end(), "features.txt") != files.end() &&
            find(files.begin(), files.end(), "rosetta_residue_scores.sc") != files.end() &&
            find(files.begin(), files.end(), "rosetta_interface_score.sc") != files.end());
}

// average the features for the ensemble
vector<double> average_features(vector<vector<double>> features) {
    vector<double> avg_features(features[0].size(), 0);
    for (size_t i = 0; i < features.size(); i++) {
        for (size_t j = 0; j < features[i].size(); j++) {
            avg_features[j] += features[i][j];
        }
    }
    for (size_t i = 0; i < avg_features.size(); i++) {
        avg_features[i] /= features.size();
    }
    return avg_features;
}

// load features for constructor
vector<EPPI::Feature> load_features(string derivative_feature_file) {
    vector<EPPI::Feature> d_features;
    ifstream f(derivative_feature_file);
    string line;
    size_t count = 0;
    while (getline(f, line)) {
        // if not in primary features list, add to the list of features
        if (find(EPPI::PrimaryFeaturesList.begin(), EPPI::PrimaryFeaturesList.end(), line) == EPPI::PrimaryFeaturesList.end()) {
            d_features.push_back(EPPI::Feature(line, count));
            count++;
        }
    }
    return d_features;
}

// Implement the method for writing to files
float PROTOCOL::EPPI_ddg(PROT::PDB* pdb, string mutation, string interface_, string output_path, bool minimize_inputs){
    // set the m_elements of the pdb for use in van der waals calculations
    pdb->set_elements();
    // the log file to keep track of calculations
    ofstream output(output_path+"/log.txt");
    // get the interface in the format required of Interface_Analyzer
    // insert a space between the interface characters break on the _ character
    string interface;
    for (size_t i = 0; i < interface_.size(); i++) {
        // if the character is _ break 
        if (interface_[i] == '_') {
            break;
        } else {
            interface += interface_[i];
            if (i != interface_.size() - 1) {
                interface += " ";
            }
        }
    }

    // make the output path
    system(("mkdir -p " + output_path).c_str());

    // gather the proteins
    vector<PROT::Protein> proteins;
    for (size_t i = 0; i < pdb->proteins(); i++) {
        proteins.push_back(*pdb->protein(i));
    }

    if (minimize_inputs) {
        ofstream min_log(output_path+"/min_log.txt");
        cout<<"Minimizing input structures"<<endl;
        // minimize the input structures
        Rosetta::Energy_Minimization(proteins, min_log, output_path);
    }

    // gather the protein pointers
    vector<PROT::Protein*> protein_ptrs;
    for (size_t i = 0; i < proteins.size(); i++) {
        protein_ptrs.push_back(&proteins[i]);
    }

    if (!already_calculated(output_path)) {
        EPPI::calculate_eppi_features(protein_ptrs, interface_, output_path, true);
    }

    // create BCProps for original pdb file
    string derived_features_file = "final_features.txt";
    vector<EPPI::Feature> derived_features = load_features(derived_features_file);
    EPPI::BCProps wt_features("wild_type", output_path, derived_features, true);
    
    // add mutation to the output path
    string mutation_output_path = output_path+"/mutation_"+ mutation;

    // make the mutation
    cout<<"Making mutation: "<<mutation<<endl;
    vector<PROT::Protein> mutated_proteins = METHODS::make_mutation(proteins, mutation, mutation_output_path);

    // run eppi features
    vector<PROT::Protein*> mutated_protein_ptrs;
    for (size_t i = 0; i < mutated_proteins.size(); i++) {
        mutated_protein_ptrs.push_back(&mutated_proteins[i]);
    }

    if (!already_calculated(mutation_output_path)) {
        EPPI::calculate_eppi_features(mutated_protein_ptrs, interface_, mutation_output_path, true);
    }
    
    // create BCProps for mutated structures
    EPPI::BCProps mutant_features(mutation, mutation_output_path, derived_features, true);

    // write the features to a file
    ofstream features_file("delta_eppi_features.txt");

    vector<double> delta_features;
    // open the model features file and write the features to the file
    string line;
    ifstream model_features(derived_features_file);
    while (getline(model_features, line)) {
        if (find(EPPI::PrimaryFeaturesList.begin(), EPPI::PrimaryFeaturesList.end(), line) != EPPI::PrimaryFeaturesList.end()) {
            size_t primary_index = find(EPPI::PrimaryFeaturesList.begin(), EPPI::PrimaryFeaturesList.end(), line) - EPPI::PrimaryFeaturesList.begin();
            double mut = mutant_features.primary_features()[primary_index];
            double wt = wt_features.primary_features()[primary_index];
            features_file<<line<<": "<<mut - wt<<endl;
        } 
        // otherwise its a derived feature
        else {
            size_t derived_index;
            for (size_t i = 0; i < derived_features.size(); i++) {
                if (derived_features[i].name() == line) {
                    derived_index = i;
                    break;
                }
            }
            double mut = mutant_features[derived_index];
            double wt = wt_features[derived_index];
            features_file<<line<<": "<<mut - wt<<endl;
        }
    }

    // run python3 PANTZ_new/models/EPPI_ddg_model.py and get the output that is printed: (ddg prediction: 0.7010752909879802)
    system("python3 PANTZ/models/EPPI_ddg_model.py --model single_state_2421 > ddg_output.txt");
    ifstream ddg_output("ddg_output.txt");
    string ddg_line;
    getline(ddg_output, ddg_line);
    // get the ddg value
    double ddg = stod(ddg_line);

    // remove the files
    system("rm ddg_output.txt");
    system("rm delta_eppi_features.txt");
    // return the ddg value
    return ddg;
    // return 0.0;
}

float PROTOCOL::EPPI_ddg_ensemble(PROT::PDB* pdb, string mutation, string interface_, string ensemble_path, bool minimize_inputs) {
    // set the m_elements of the pdb for use in van der waals calculations
    pdb->set_elements();
    cout<<"Calculating features for ensemble mutations: "<<mutation<<endl;
    string derived_features_file = "final_features.txt";

    // calculate the EPPU features for the original pdb
    vector<PROT::Protein*> proteins;
    for (size_t i = 0; i < pdb->proteins(); i++) {
        proteins.push_back(pdb->protein(i));
    }
    if (!already_calculated(ensemble_path)) {
        EPPI::calculate_eppi_features(proteins, interface_, ensemble_path, true);
    }
    // create BCProps for original pdb file
    vector<EPPI::Feature> derived_features = load_features(derived_features_file);
    EPPI::BCProps wt_features("wild_type", ensemble_path, derived_features, true);

    // calculate the average features for the ensemble
    vector<string> dir = METHODS::listdir(ensemble_path+"/");
    // keep only the pdb files
    vector<string> ensemble_files;
    for (size_t i = 0; i < dir.size(); i++) {
        if (dir[i].find(".pdb") != string::npos) {
            ensemble_files.push_back(dir[i]);
        }
    }
    vector<vector<double>> avg_features;
    // calculate the features for each pdb file
    for (size_t i = 0; i < ensemble_files.size(); i++) {
        vector<PROT::Protein> proteins;
        // load the pdb file
        PROT::PDB pdb(ensemble_path+"/"+ensemble_files[i]);
        for (size_t j = 0; j < pdb.proteins(); j++) {
            proteins.push_back(*pdb.protein(j));
        }
        // remove .pdb
        ensemble_files[i] = ensemble_files[i].substr(0, ensemble_files[i].size()-4);
        if (minimize_inputs) {
            ofstream min_log(ensemble_path+"/"+ensemble_files[i]+"/min_log.txt");
            cout<<"Minimizing input structures"<<endl;
            // minimize the input structures
            Rosetta::Energy_Minimization(proteins, min_log, ensemble_path+"/"+ensemble_files[i]);
        }
        // gather the protein pointers
        vector<PROT::Protein*> protein_ptrs;
        for (size_t j = 0; j < proteins.size(); j++) {
            protein_ptrs.push_back(&proteins[j]);
        }
        if (!already_calculated(ensemble_path+"/"+ensemble_files[i])) {
            string output_path = ensemble_path+"/"+ensemble_files[i];
            EPPI::calculate_eppi_features(protein_ptrs, interface_, output_path, true);
        }
        // create BCProps for the mutated structures
        EPPI::BCProps features(ensemble_files[i], ensemble_path+"/"+ensemble_files[i], true);
        avg_features.push_back(features.base_features());
    }

    // average the base features of the ensemble
    vector<double> avg = average_features(avg_features);

    vector<string> ensemble_base_features;
    ensemble_base_features.push_back(ensemble_path);
    for (size_t i = 0; i < avg.size(); i++) {
        ensemble_base_features.push_back(to_string(avg[i]));
    }

    // create an ensemble bcprops
    EPPI::BCProps ensemble_features(ensemble_path, ensemble_base_features, derived_features);

    // write the features to a file
    ofstream features_file("delta_eppi_features.txt");

    vector<double> delta_features;
    // open the model features file and write the features to the file
    string line;
    ifstream model_features(derived_features_file);
    while (getline(model_features, line)) {
        if (find(EPPI::PrimaryFeaturesList.begin(), EPPI::PrimaryFeaturesList.end(), line) != EPPI::PrimaryFeaturesList.end()) {
            size_t primary_index = find(EPPI::PrimaryFeaturesList.begin(), EPPI::PrimaryFeaturesList.end(), line) - EPPI::PrimaryFeaturesList.begin();
            double mut = ensemble_features.primary_features()[primary_index];
            double wt = wt_features.primary_features()[primary_index];
            features_file<<line<<": "<<mut - wt<<endl;
        } 
        // otherwise its a derived feature
        else {
            size_t derived_index;
            for (size_t i = 0; i < derived_features.size(); i++) {
                if (derived_features[i].name() == line) {
                    derived_index = i;
                    break;
                }
            }
            double mut = ensemble_features[derived_index];
            double wt = wt_features[derived_index];
            features_file<<line<<": "<<mut - wt<<endl;
        }
    }

    // run python3 PANTZ_new/models/EPPI_ddg_model.py and get the output that is printed: (ddg prediction: 0.7010752909879802)
    system("python3 PANTZ/models/EPPI_ddg_model.py --model ensemble_2421 > ddg_output.txt");
    ifstream ddg_output("ddg_output.txt");
    string ddg_line;
    getline(ddg_output, ddg_line);
    // get the ddg value
    double ddg = stod(ddg_line);    
    
    // remove the files
    system("rm ddg_output.txt");
    system("rm delta_eppi_features.txt");
    // return the ddg value
    return ddg;
}

