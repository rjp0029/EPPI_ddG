/* Created by clay at Auburn University.
*
* This file implements the functions to make mutation to a reisdue
* This is based on work done in the paper by
* Richard et al. 2025.
*
*/

// This file is supposed to be loaded by Methods.h
#ifndef Methods_Loading_Status
#error Muation Methods must be included by Methods.h
#endif

// mutation function (mutation -> DA26A <old><chain><residue number><new>) 
PROT::PDB METHODS::make_mutation(PROT::PDB * pdb, string& mutation, string& output_path){
    char mutant_res = mutation[0];
    char mutant_chain = mutation[1];
    int mutant_res_num = std::stoi(mutation.substr(2, mutation.size() - 3));
    char new_res_char = mutation.back();

    // make the output path
    system(("mkdir -p " + output_path).c_str());
 
    // convert new res to 3 letter code
    string new_res;
    for (size_t i = 0; i < PROT::AA1.size(); i++) {
        if (PROT::AA1[i][0] == new_res_char) {
            new_res = PROT::AA3[i];
            break;
        }
    }

    ofstream output (output_path+"/log.txt");
    vector<PROT::Protein> proteins;
    for (size_t i = 0; i < pdb->proteins(); i++) {
        proteins.push_back(*pdb->protein(i));
    }
    // make the mutation
    PROT::Residue* res;
    for (size_t i = 0; i < proteins.size(); i++) {
        for (size_t j = 0; j < proteins[i].size(); j++) {
            if (proteins[i].operator()(j, ' ', true)->protein() == mutant_chain && proteins[i].operator()(j, ' ', true)->AA1() == mutant_res && proteins[i].operator()(j, ' ', true)->number() == mutant_res_num) {
                res = proteins[i].operator()(j, ' ', true);
                break;
            }
        }
    }
    // rename the residue
    res->remove_sidechain();
    res->rename(new_res);

    // run a minimization with all residues fixed except the mutated residue
    vector<int> fixed_residues;
    int fix_res = Rosetta::get_pose_numbering(pdb, mutant_chain, mutant_res_num);
    fixed_residues.push_back(fix_res);
    Rosetta::Energy_Minimization_fixed_res(proteins, output, output_path, fixed_residues);

    ofstream mutated_file_minimized (output_path+"/mutation_" + mutation + "_minimized.pdb");
    for (size_t i = 0; i < proteins.size(); i++) {
        mutated_file_minimized<<proteins[i].str()<<endl;
    }
    mutated_file_minimized.close();
    return PROT::PDB(output_path+"/mutation_" + mutation + "_minimized.pdb");
}

// mutation function (mutation -> DA26A <old><chain><residue number><new>) (interface -> A_BC <chain(s)>_<chain(s)>)
vector<PROT::Protein> METHODS::make_mutation(vector<PROT::Protein>& proteins, string& mutation, string& output_path){
    char mutant_res = mutation[0];
    char mutant_chain = mutation[1];
    int mutant_res_num = std::stoi(mutation.substr(2, mutation.size() - 3));
    char new_res_char = mutation.back();

    // make the output path
    system(("mkdir -p " + output_path).c_str());

    // convert new res to 3 letter code
    string new_res;
    for (size_t i = 0; i < PROT::AA1.size(); i++) {
        if (PROT::AA1[i][0] == new_res_char) {
            new_res = PROT::AA3[i];
            break;
        }
    }

    ofstream output (output_path+"/log.txt");
    // make the mutation
    PROT::Residue* res;
    for (size_t i = 0; i < proteins.size(); i++) {
        for (size_t j = 0; j < proteins[i].size(); j++) {
            if (proteins[i].operator()(j, ' ', true)->protein() == mutant_chain && proteins[i].operator()(j, ' ', true)->AA1() == mutant_res && proteins[i].operator()(j, ' ', true)->number() == mutant_res_num) {
                res = proteins[i].operator()(j, ' ', true);
                break;
            }
        }
    }
    // rename the residue
    res->remove_sidechain();
    res->rename(new_res);

    // run a minimization with all residues fixed except the mutated residue
    vector<int> fixed_residues;
    int fix_res = Rosetta::get_pose_numbering(proteins, mutant_chain, mutant_res_num);
    fixed_residues.push_back(fix_res);
    Rosetta::Energy_Minimization_fixed_res(proteins, output, output_path, fixed_residues);
    // write to a file
    ofstream mutated_file_minimized (output_path+"/mutation_" + mutation + "_minimized.pdb");
    for (size_t i = 0; i < proteins.size(); i++) {
        mutated_file_minimized<<proteins[i].str()<<endl;
    }
    mutated_file_minimized.close();
    return proteins;
}