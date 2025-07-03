/* created by the PROTEINPANT(z) lab at Auburn University 2025
    This is the file that implements the EPPI::Residue class
*/

// preproccessor directives to make sure the files is included by the
// correct file
#ifndef EPPI_Loading_Status
#error Residue.h must be included by BCProps.h
#endif


// these classes are separate from the Pantz2 
class EPPI::Residue {
    public:
        Residue(string name, int number, char chain, 
            int pre_free_rotamers, int bound_free_rotamers, string prestable, string boundstable);
        string name;
        int number;
        char chain;
        int pre_free_rotamers = 0;
        int bound_free_rotamers = 0;
        bool prestable;
        bool boundstable;
};

// this is the constructor for the EPPIResidue class
EPPI::Residue::Residue(string name, int number, char chain, 
        int pre_free_rotamers, int bound_free_rotamers, string prestable, string boundstable) {
    this->name = name;
    this->number = number;
    this->chain = chain;
    this->pre_free_rotamers = pre_free_rotamers;
    this->bound_free_rotamers = bound_free_rotamers;
    if (prestable == "prestable") {
        this->prestable = true;
    } else {
        this->prestable = false;
    }
    if (boundstable == "bound_stable") {
        this->boundstable = true;
    } else {
        this->boundstable = false;
    }
}