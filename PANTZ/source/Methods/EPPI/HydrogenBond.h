/* created by the PROTEINPANT(z) lab at Auburn University 2025
    This is the file that implements the EPPI::HydrogenBond class
*/

// preproccessor directives to make sure the files is included by the
// correct file
#ifndef EPPI_Loading_Status
#error HydrogenBond.h must be included by EPPI.h
#endif

// this is the hydrogen bond class that inherits from the EPPIInteraction class
class EPPI::HydrogenBond : public EPPI::Interaction {
    public:
        HydrogenBond(string line);
        double distance;
        double dha_angle;
        double daa_angle;
        bool is_real() {return this->residue1->boundstable && this->residue2->boundstable;};
        string str() {
            return "Hydrogen Bond "+ residue1->name + to_string(residue1->number) + residue1->chain + " " +
                to_string(residue1->pre_free_rotamers) + " " + to_string(residue1->bound_free_rotamers) +
                " - " + residue2->name + to_string(residue2->number) + residue2->chain + " " +
                to_string(residue2->pre_free_rotamers) + " " + to_string(residue2->bound_free_rotamers) +
                " " + to_string(distance) + " " + to_string(dha_angle) + " " + to_string(daa_angle) +
                " " + to_string(energy);
        }
};

// this is the constructor for the EPPIHydrogenBond class
EPPI::HydrogenBond::HydrogenBond(string line) : EPPI::Interaction(line) {
    vector<string> words;
    Text::split(words, line);
    this->distance = stof(words[16]);
    this->dha_angle = stof(words[17]);
    this->daa_angle = stof(words[18]);
}