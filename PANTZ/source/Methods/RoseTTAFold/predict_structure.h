/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements energy minimization functions for RoseTTAFold */

// Make sure this file is being included from the RoseTTAFold.h header file
#ifndef RoseTTAFold_Loading_Status
#error RoseTTAFold functions must be included from RoseTTAFold.h
#endif

// The function that actually runs a RoseTTAFold structure prediction
PROT::Protein * RoseTTAFold::predict_structure(string & primary_sequence, char & protein_name) {
    // create a fasta file with the primary sequence for RoseTTAFold
    ofstream fasta_file("prediction.fa");
    fasta_file<<primary_sequence;
    fasta_file.close();
    // // print the details of the prediction
    cout<<"predicting the structure of the PROT::Protein with the primary sequence: \n"<<primary_sequence<<endl;
    cout<<"\nThis may take a while..."<<endl;
    // create the output name using the string of the protein name
    string rf_output_name = "./rf_output_" + string(1, protein_name);
    // the command to run the pyrosetta mode of RoseTTAFold
    string command = string(ROSETTAFOLD_PATH) + "/run_pyrosetta_ver.sh ./prediction.fa " + rf_output_name;
    system(command.c_str());
    // get the protein that was predicted
    PROT::Protein predicted_protein("model_1.crderr.pdb", "./rf_output_" + string(1, protein_name) + "/model/");
    // rename the predicted protein to the name of the original
    predicted_protein.set_name(protein_name);
    // return the pointer to the predicted protein
    return new PROT::Protein(predicted_protein);
}
