/* Created by clay at Auburn University.
*
* This file implements the protein aligning method using Singular Value
* Decomposition. 
*
*/

// This file is supposed to be loaded by Methods.h
#ifndef Methods_Loading_Status
#error General Methods must be included by Methods.h
#endif

// a function to get the atoms of the specified type
vector<PROT::Atom*> get_atoms (vector<PROT::Residue*> residues, const string& how) {
    vector<PROT::Atom*> atoms;
    if (how == "all") {
        for (auto residue : residues) {
            for (size_t j = 0; j < residue->size(); j++) {
                atoms.push_back(residue->get_atom(j));
            }
        }
    } else if (how == "backbone") {
        for (auto residue : residues) {
            atoms.push_back(residue->get_atom("N"));
            atoms.push_back(residue->get_atom("CA"));
            atoms.push_back(residue->get_atom("C"));
            atoms.push_back(residue->get_atom("O"));
        }
    } else if (how == "CA") {
        for (auto residue : residues) {
            atoms.push_back(residue->get_atom("CA"));
        }
    } else {
        string error = "The atom type " + how + " is not recognized.";
        throw PANTZ_error(error);
    }
    return atoms;
}

// get common residues from the two proteins
vector<vector<PROT::Residue*>> get_common_residues (PROT::Protein* reference, PROT::Protein* target) {
    vector<PROT::Residue*> reference_residues;
    vector<PROT::Residue*> target_residues;
    for (size_t i = 0; i < reference->size(); i++) {
        // if the residue is missing, print it and skip to the next residue
        if (!reference->operator()(i, ' ', true)->is_present()) {
            cout<<"Missing residue: "<<reference->operator()(i, ' ', true)->name()<<endl;
            continue;
        }
        for (size_t j = 0; j < target->size(); j++) {
            if (!target->operator()(j, ' ', true)->is_present()) {
                cout<<"Missing residue: "<<target->operator()(j, ' ', true)->name()<<endl;
                continue;
            }
            if ((reference->operator()(i, ' ', true)->name() == target->operator()(j, ' ', true)->name()) and
                (reference->operator()(i, ' ', true)->internal_number() == target->operator()(j, ' ', true)->internal_number())) {
                // (reference->operator()(i, ' ', true)->protein() == target->operator()(j, ' ', true)->protein())) {
                reference_residues.push_back(reference->operator()(i, ' ', true));
                target_residues.push_back(target->operator()(j, ' ', true));
                // cout<<"Common residue found: "<<reference->operator()(i, ' ', true)->name()<<endl;
            }
        }
    }
    vector<vector<PROT::Residue*>> common_residues;
    common_residues.push_back(reference_residues);
    common_residues.push_back(target_residues);
    // if no common residues are found, throw an error
    if (common_residues[0].size() == 0) {
        string error = "No common residues found between the two proteins. Must have same name, number, and protein.";
        throw PANTZ_error(error);
    }
    return common_residues;
}

// the align method that aligns the target onto the reference using rotation matrix 
// calculated from by SVD of the covariance matrix of the two sets of coordinates
// aligning the second set of coordinates onto the first set of coordinates
void METHODS::align (PROT::Protein* reference, PROT::Protein* target, vector<vector<PROT::Atom*>> atoms) {
    // get the coordinates of the reference and target atoms
    PROT::Matrix ref_coords = PROT::Matrix(atoms[0]);
    PROT::Matrix target_coords = PROT::Matrix(atoms[1]);
    // align the target matrix onto the reference matrix 
    PROT::Matrix rot_matrix = ref_coords.rotation_matrix(&target_coords);
    // get the centroids of the reference and target proteins
    PROT::Matrix ref_centroid = ref_coords.average(true);
    PROT::Matrix target_centroid = target_coords.average(true);
    // apply the movement to all atoms in the target protein
    target->move(target_centroid, '-');
    target->rotate(rot_matrix);
    target->move(ref_centroid, '+');
    return;
}


// Implement the function to align the target protein onto the reference protein using the 
// specified atoms (i.e., all, backbone, or CA
void METHODS::global_align (PROT::Protein* reference, PROT::Protein* target, const string& how = "all") {
    // check that the sizes are equal, print the sizes if not 
    if (reference->size() != target->size()) {
        cout<<"Reference protein size: "<<reference->size()<<endl;
        cout<<"Target protein size: "<<target->size()<<endl;
        string error = "The two proteins must have the same number of residues.";
        throw PANTZ_error(error);
    }
    cout<<"Protein sizes are equal at "<<reference->size()<<endl;
    // get common residues
    vector<vector<PROT::Residue*>> common_residues = get_common_residues(reference, target);
    cout<<"Common residues found: "<<common_residues[0].size()<<endl;
    // get the atoms of the specified type 
    vector<PROT::Atom*> ref_atoms;
    vector<PROT::Atom*> target_atoms;
    ref_atoms = get_atoms(common_residues[0], how);
    target_atoms = get_atoms(common_residues[1], how);
    // align the target onto the reference using the atoms
    align(reference, target, {ref_atoms, target_atoms});
    return;
}