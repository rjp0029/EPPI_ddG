/* Created by clay at Auburn University.
 *
 * This file implements function to get pose numbering from Rosetta */

// Make sure this file is being included from the Rosetta.h header file
#ifndef Rosetta_Loading_Status
#error Rosetta functions must be included from Rosetta.h
#endif

// method to get the pose numbering
int Rosetta::get_pose_numbering(string pdb_file, char chain, int pdb_res_num){
    size_t pose_num = 0;
    // Rosetta numbering of residues in a protein structure, sometimes called Pose numbering, is distinct from the residue numbering in the input PDB files ("PDB numbering"). 
    // Rosetta numbering always begins at 1 for the first residue and increases by one for each residue, ignoring chain designation.
    PROT::PDB pdb(pdb_file);
    for (size_t i = 0; i < pdb.proteins(); i++) {
        for (size_t j = 0; j < pdb.protein(i)->size(); j++) {
            pose_num += 1;
            if (pdb.protein(i)->operator()(j, ' ', true)->protein() == chain && pdb.protein(i)->operator()(j, ' ', true)->number() == pdb_res_num) {
                break;
            }
        }
    }
    return pose_num;
}

// method to get the pose numbering
int Rosetta::get_pose_numbering(PROT::PDB* pdb, char chain, int pdb_res_num){
    size_t pose_num = 0;
    // Rosetta numbering of residues in a protein structure, sometimes called Pose numbering, is distinct from the residue numbering in the input PDB files ("PDB numbering"). 
    // Rosetta numbering always begins at 1 for the first residue and increases by one for each residue, ignoring chain designation.
    for (size_t i = 0; i < pdb->proteins(); i++) {
        for (size_t j = 0; j < pdb->protein(i)->size(); j++) {
            pose_num += 1;
            if (pdb->protein(i)->operator()(j, ' ', true)->protein() == chain && pdb->protein(i)->operator()(j, ' ', true)->number() == pdb_res_num) {
                break;
            }
        }
    }
    return pose_num;
}

// method to get the pose numbering
int Rosetta::get_pose_numbering(vector<PROT::Protein>& proteins, char chain, int pdb_res_num){
    size_t pose_num = 0;
    // Rosetta numbering of residues in a protein structure, sometimes called Pose numbering, is distinct from the residue numbering in the input PDB files ("PDB numbering"). 
    // Rosetta numbering always begins at 1 for the first residue and increases by one for each residue, ignoring chain designation.
    for (size_t i = 0; i < proteins.size(); i++) {
        for (size_t j = 0; j < proteins[i].size(); j++) {
            pose_num += 1;
            if (proteins[i].operator()(j, ' ', true)->protein() == chain && proteins[i].operator()(j, ' ', true)->number() == pdb_res_num) {
                break;
            }
        }
    }
    return pose_num;
}
