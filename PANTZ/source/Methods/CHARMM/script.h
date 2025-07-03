/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements methods that make parts of CHARMM scripts. */

// This file must be included from CHARMM.h
#ifndef CHARMM_Loading_Status
#error CHARMM methods must be included by CHARMM.h
#endif

// Text that sets warning and bomb levels
void CHARMM::warning_bomb (string& script) {
    script += "wrnl -2\nbomb -2\n\n";
}

// Text that loads topology and parameter files
void CHARMM::load_topology_parameter (string& script) {
    script += "! Load the Topology File\n"
              "open read unit 10 form name "
              "/home/shared/rjp0029_lab/charmm/toppar/top_all36_prot.rtf"
              " card\nread rtf card unit 10\nclose unit 10\n\n"
              "! Load the Parameter File\n"
              "open read unit 10 form name "
              "/home/shared/rjp0029_lab/charmm/toppar/par_all36_prot.prm"
              " card\nread para card unit 10\nclose unit 10\n\n";
}

// Text that creates missing atoms
void CHARMM::add_missing_atoms (string& script) {
    script += "! Add missing atoms\n"
              "ic fill preserve\n"
              "ic param\n"
              "ic build\n"
              "hbuild\n\n";
}

// Text that runs an energy minimization
void CHARMM::energy_minimization (string& script, const bool harmonic,
                                  const bool fixedBackbone) {
    // Update the comment in the script appropriately
    if ((!harmonic) && (!fixedBackbone)) {
        script += "! All Atom Energy Minimization\n";}
    else if ((harmonic) && (!fixedBackbone)) {
        script += "! Harmonic Backbone Energy Minimization\n";}
    else if ((!harmonic) && (fixedBackbone)) {
        script += "! Fixed Backbone Energy Minimization\n";}
    else {
        string error = "Algorithm Error: Invalid inputs to "
                       "CHARMM::energy_minimization\n";
        throw PANTZ_error (error);}
    // If a harmonic constraint is being used
    if (harmonic) {
        script += "cons harm force 10 sele ( type CA .or. type C .or. type N) end\n";}
    // if a fixed backbone constraint is being used
    else if (fixedBackbone) {
        script += "cons fix sele ( type CA .or. type C .or. type N) end\n";}
    // Update the rest of the energy minimization code
    script += "nbon nbxm 5\n"
              "skip all excl angl bond dihe impr urey elec vdw";
    if (harmonic) {script += " harm";}
    script += "\nmini abnr nstep 5000 nprint 10 -\n"
              "tolgrd 0.01 tolenr 0.0001 tolstp 0.00\n\n";
}

// Text that tells charmm to output proteins to a file
void CHARMM::output_proteins (string& script, vector<PROT::Protein>& prots) {
    // Do this for each protein
    for(size_t i=0; i<prots.size(); ++i) {
        // Get the name of the file to write the protein to
        string fileName = make_protein_name (prots[i], false);
        // Add a comment to the script
        script += "! Output Protein ";
        script += prots[i].name();
        script += "\nopen write unit 10 name " + fileName + " card\n"
                  "write coor sele segi pr";
        string name = "";
        name += prots[i].name();
        Text::lower(name);
        script += name + " end pdb unit 10 card\n"
                  "close unit 10\n\n";}
}
