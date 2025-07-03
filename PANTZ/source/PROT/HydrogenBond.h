/* Created by the clay at Auburn University.
 *
 * This file contains the declaration of the HydrogenBond struct. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_HydrogenBond_Guard
#define Proteins_HydrogenBond_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error HydrogenBond.h must be included by Proteins.h
#endif

// confirm that the Residue header files have been included
#ifndef Proteins_Residue_Guard
#error HydrogenBond.h must be included after Residue.h
#endif

// define the hydrogen bond struct
struct PROT::HydrogenBond{
    PROT::Residue* donor_residue;
    bool donor_backbone = false;
    size_t donor_free_rot_pre;
    size_t donor_free_rot_bound;
    PROT::Residue* acceptor_residue;
    bool acceptor_backbone = false;
    size_t acceptor_free_rot_pre;
    size_t acceptor_free_rot_bound;
    PROT::Atom* donor_atom;
    PROT::Atom* hydrogen;
    PROT::Atom* acceptor_atom;
    PROT::Atom* antecedent;
    float distance;
    float DHA_angle;
    float DAAn_angle;
    // equal operator that compares the donor and acceptor residues and backbones
    bool operator==(const HydrogenBond& other){
        return (donor_residue->name() == other.donor_residue->name() && 
                acceptor_residue->name() == other.acceptor_residue->name() &&
                donor_backbone == other.donor_backbone && 
                acceptor_backbone == other.acceptor_backbone &&
                donor_atom->name() == other.donor_atom->name() && 
                acceptor_atom->name() == other.acceptor_atom->name() &&
                hydrogen->name() == other.hydrogen->name());
    }
};

// End the header guard from the start of the file
#endif
