/* Created by clay at Auburn University.
 *
 * This file implements the protocol of filling in missing residues 
   in a crystal structure. It does so by predicting the structure of 
   the protein object using rosettafold, then aligning the structures
   with the Kabsch algorithm (method) and then splicing in the missing
   residues from the predicted structure. Finally the resulting structure
   is energy minimized to correct for sterics and backbone anomalies. */

// This file is supposed to be included by Protocol.h
#ifndef Protocol_Loading_Status
#error General Protocol functions must be included by Protocol.h
#endif

// get the gap residues (first and last residue of the gap, inclusive of the end points) (there can be multiple gaps)
vector<vector<PROT::Residue*>> get_gaps (PROT::Protein* protein) {
    vector<vector<PROT::Residue*>> gap_residues;
    vector<PROT::Residue*> gap;
    // go through residues in the protein
    for (size_t i = 0; i < protein->size(); i++) {
        // if residue is not present
        if (!protein->operator()(i, ' ', true)->is_present()) {
            // add it to the gap
            gap.push_back(protein->operator()(i, ' ', true));
        } else {
            // if the gap is not empty, add it to the list of gaps
            if (gap.size() > 0){
                gap_residues.push_back(gap);
                gap.clear();
            }
        }
    }
    // get the edge case where there is a terminal gap
    if (!gap.empty()) {
        gap_residues.push_back(gap);
    }
    return gap_residues;
}

// fill the gaps in a protein object if the predicted structure is already known
void PROTOCOL::align_and_splice_residues(PROT::Protein * protein, PROT::Protein * predicted_protein, bool terminal) {
    // // check for non amino acid residues
    // for (size_t j = 0; j < protein->size(); j++){
    //     if (!protein->operator()(j, ' ', true)->is_amino_acid()){
    //         cout<<"Non amino acid residue found: "<<protein->operator()(j, ' ', true)->name()<<endl;
    //         return;
    //     }
    // }
    // rename the predicted protien to match the original protein
    predicted_protein->set_name(protein->name());
    // get the gap residues in the original protein
    vector<vector<PROT::Residue*>> gap_residues = get_gaps(protein);
    PROT::Residue predicted_residue;
    // go through each gap and align the predicted protein to the original protein based 
    // on the atoms of the gap residues (backbone only) and then replace the gap residues 
    // with the corresponding residues in the predicted protein
    for (size_t i = 0; i < gap_residues.size(); i++){
        // get the atoms of the residues before and after the gap
        cout<<"Gap "<<i+1<<": "<<gap_residues[i][0]->internal_number()<<" - "<<gap_residues[i][gap_residues[i].size()-1]->internal_number()<<" identified"<<endl;
        // if terminal is false, then skip gaps that are contnious with the first or last residue
        if (!terminal){
            if (gap_residues[i][0]->internal_number() == 1 || gap_residues[i][gap_residues[i].size()-1]->internal_number() == protein->size()){
                cout<<"Skipping terminal gap, enable terminal=true to fill terminal gaps"<<endl;
                continue;
            }
        }
        // create vectors of atoms for the reference and target proteins for alinging the anchor residues
        vector<PROT::Atom*> reference_atoms;
        vector<PROT::Atom*> target_atoms;
        // get the anchor residue backbone atoms to attach to
        int start = gap_residues[i][0]->internal_number()-1;
        int end = gap_residues[i][gap_residues[i].size()-1]->internal_number()+1;
        if (start > 0){
            protein->operator()(start-1, ' ', true)->select_atoms(reference_atoms, "backbone heavy");
            predicted_protein->operator()(start-1, ' ', true)->select_atoms(target_atoms, "backbone heavy");
        }
        if (end < protein->size()){
            protein->operator()(end-1, ' ', true)->select_atoms(reference_atoms, "backbone heavy");
            predicted_protein->operator()(end-1, ' ', true)->select_atoms(target_atoms, "backbone heavy");
        }
        // align the target protein to the reference protein
        METHODS::align(protein, predicted_protein, {reference_atoms, target_atoms});
        // go through the reference protein (protein_nonhetatm):
        // if the residue is missing, then replace it with the 
        // // residue in the target protein with the corresponding internal number
        for (size_t j = 0; j < gap_residues[i].size(); j++){
            int original_number = gap_residues[i][j]->number();
            // replace the residue
            predicted_residue = predicted_protein->operator()(gap_residues[i][j]->internal_number() - 1, ' ', true);
            // cout<<"predicted residue: "<<predicted_residue.name()<<" "<<predicted_residue.number()<<endl;
            predicted_residue.set_number(original_number, ' ', false);
            // cout<<"renumbered to: "<<predicted_residue.name()<<" "<<predicted_residue.number()<<endl;
            // cout<<"Changed from "<<gap_residues[i][j]->name()<<" "<<gap_residues[i][j]->number();
            // *gap_residues[i][j] = *predicted_residue;
            // destruct the original residue from the gap residues and replace it with the predicted residue
            gap_residues[i][j]->~Residue();
            new (gap_residues[i][j]) PROT::Residue(predicted_residue);
            // cout << " to " << gap_residues[i][j]->name() << " " << gap_residues[i][j]->number() << endl;
        }
    }
    protein->renumber_atoms(1);
    cout<<"Done filling gaps in "<<protein->name()<<endl;
}

// Fill the gaps in a Protein Object predicting the structure along the way
void PROTOCOL::fill_gaps(PROT::Protein * protein) {
    // get the name of the protein
    char protein_name = protein->name();
    // create a new protein using a vector of residues and populate with only the hetatm lines
    vector<PROT::Residue> residues;
    for (size_t j = 0; j < protein->size(); j++){
        // if residue is empty, add it to the list of residues
        if (protein->operator()(j, ' ', true)->size() == 0){
            residues.push_back(*protein->operator()(j, ' ', true));
            continue;
        }
        // if the residue does not contain a HETATM line, add it to the list of residues
        if (protein->operator()(j, ' ', true)->get_atom(0)->type() != "HETATM"){
            residues.push_back(*protein->operator()(j, ' ', true));
        }
    }
    // create a new protein object with the non HETATM residues
    PROT::Protein protein_nonhetatm(residues);
    // get the primary sequence of the protein
    string sequence = protein_nonhetatm.fasta();
    // predict the protein structure
    PROT::Protein * predicted_protein = RoseTTAFold::predict_structure(sequence, protein_name);
    // get the gaps and fill with the predicted structure
    align_and_splice_residues(protein, predicted_protein);
}
