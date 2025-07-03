/* created by the PROTEINPANT(z) lab at Auburn University 2025
    This is the file that implements the EPPI::SaltBridge class
*/

// preproccessor directives to make sure the  files is included by the
// correct file
#ifndef EPPI_Loading_Status
#error SaltBridge.h must be included by BCProps.h
#endif



// this is the class for a salt bridge that inherits from the EPPIInteraction class
class EPPI::SaltBridge : public EPPI::Interaction {
    public:
        SaltBridge(string line);
        double distance;
        bool is_real() {return true;};
        string str() {
            return "Salt Bridge " + residue1->name + to_string(residue1->number) + residue1->chain + " " +
                to_string(residue1->pre_free_rotamers) + " " + to_string(residue1->bound_free_rotamers) +
                " - " + residue2->name + to_string(residue2->number) + residue2->chain + " " +
                to_string(residue2->pre_free_rotamers) + " " + to_string(residue2->bound_free_rotamers) +
                " " + to_string(distance) + " " + to_string(energy);
        }
};

// this is the constructor for the EPPISaltBridge class
EPPI::SaltBridge::SaltBridge(string line) : EPPI::Interaction(line) {
    vector<string> words;
    Text::split(words, line);
    this->distance = stof(words[16]);
}