/* created by the PROTEINPANT(z) lab at Auburn University 2025
    This is the file that implements the EPPI::Hydrophobic class
*/

// preproccessor directives to make sure the file is included by the
// correct file
#ifndef EPPI_Loading_Status
#error Hydrophobic.h must be included by BCProps.h
#endif


// this is the class for a hydrophobic interaction that inherits from the EPPIInteraction class
class EPPI::Hydrophobic : public Interaction {
    public:
        Hydrophobic(string line);
        double bsasa;
        bool is_real() {return true;};
        string str() {
            return "Hydrophobic Interaction " + residue1->name + to_string(residue1->number) + residue1->chain + " " +
                to_string(residue1->pre_free_rotamers) + " " + to_string(residue1->bound_free_rotamers) +
                " - " + residue2->name + to_string(residue2->number) + residue2->chain + " " +
                to_string(residue2->pre_free_rotamers) + " " + to_string(residue2->bound_free_rotamers) +
                " " + to_string(bsasa) + " " + to_string(energy);
        }
};

// this is the constructor for the EPPIHydrophobic class
EPPI::Hydrophobic::Hydrophobic(string line) : EPPI::Interaction(line) {
    vector<string> words;
    Text::split(words, line);
    this->bsasa = stof(words[16]);
}