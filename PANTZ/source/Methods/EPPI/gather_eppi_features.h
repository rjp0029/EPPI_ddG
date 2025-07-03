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

void EPPI::gather_eppi_features(string path, vector<string> dir_names, 
    string derivative_feature_file, string csv_file, bool verbose){

    // remove the csv
    system(("rm " + csv_file).c_str());

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

    ofstream f_csv(csv_file);
    // write the features to a csv file
    f_csv<<"Name,";
    for (size_t i=0; i<EPPI::PrimaryFeaturesList.size(); ++i) {
        f_csv<<EPPI::PrimaryFeaturesList[i]<<",";
    }
    for (size_t i=0; i<d_features.size(); ++i) {
        f_csv<<d_features[i].name()<<",";
    }
    f_csv<<endl;

    string output_dir;
    for (size_t i=0; i<dir_names.size(); ++i) {
        output_dir = path + dir_names[i];
        if (verbose) {
            cout<<"Gathering features for "<<dir_names[i]<<endl;
        }
        try {
            EPPI::BCProps features(dir_names[i], output_dir, d_features, verbose);
            features.append_to_csv(csv_file);
        } // if a PANTZ_error is thrown, print the error message but continue to the next complex
        catch (PANTZ_error& e) {
            cout<<e.what()<<endl;
        }
        // EPPI::BCProps features(dir_name[i], output_dir, d_features, verbose);
        // // cout<<"constructed BCProps"<<endl;
        // features.append_to_csv(csv_file);s
    }
}
