/* created by the PROTEINPANT(z) lab at Auburn University 2025
    This is the file that implements the EPPI::Interaction class
*/

// preproccessor directives to make sure the files is included by the
// correct file
#ifndef EPPI_Loading_Status
#error Interaction.h must be included by BCProps.h
#endif


// this is the class for an EPPI interaction
class EPPI::Interaction {
    public:
        Interaction(string line);
        double energy = 0.0;
        EPPI::Residue* residue1;
        EPPI::Residue* residue2;
        void set_energy(map<string, double>& energy_map);
        // the string method is used to print the interaction
        string str() {
            return residue1->name + to_string(residue1->number) + residue1->chain + " " +
                to_string(residue1->pre_free_rotamers) + " " + to_string(residue1->bound_free_rotamers) +
                " - " + residue2->name + to_string(residue2->number) + residue2->chain + " " +
                to_string(residue2->pre_free_rotamers) + " " + to_string(residue2->bound_free_rotamers);
        }
};

// this is the constructor for the EPPIInteraction class
EPPI::Interaction::Interaction(string line) {
    vector<string> words;
    Text::split(words, line);
    this->residue1 = new EPPI::Residue(words[2].substr(0, 3), stoi(words[2].substr(4, words[2].size() - 6)), 
        words[2][words[2].size() - 1], stoi(words[6]), 
        stoi(words[7]), words[3], words[4]);
    this->residue2 = new EPPI::Residue(words[9].substr(0, 3), stoi(words[9].substr(4, words[9].size() - 6)), 
        words[9][words[9].size() - 1], stoi(words[13]), 
        stoi(words[14]), words[10], words[11]);
}

// set the energy of the interaction from the energy map
void EPPI::Interaction::set_energy(map<string, double>& energy_map) {
    // the key could be 1->2 or 2->1
    string key1 = to_string(residue1->number) + residue1->chain + residue1->name + " " + 
        to_string(residue2->number) + residue2->chain + residue2->name;
    string key2 = to_string(residue2->number) + residue2->chain + residue2->name + " " +
        to_string(residue1->number) + residue1->chain + residue1->name;
    if (energy_map.find(key1) != energy_map.end()) {
        this->energy = energy_map[key1];
    } else if (energy_map.find(key2) != energy_map.end()) {
        this->energy = energy_map[key2];
    } else {
        throw PANTZ_error("Error in set_energy: key not found for: " + key1 + " or " + key2);
    }
}


