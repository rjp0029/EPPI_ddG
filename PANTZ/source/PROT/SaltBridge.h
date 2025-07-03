/* Created by the clay at Auburn University.
 *
 * This file contains the declaration of the SaltBridge struct. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_SaltBridge_Guard
#define Proteins_SaltBridge_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error SaltBridge.h must be included by Proteins.h
#endif

// Define the SaltBridge struct
struct PROT::SaltBridge{
    PROT::Residue* donor_residue;
    bool donor_backbone = false;
    size_t donor_free_rot_pre;
    size_t donor_free_rot_bound;
    PROT::Residue* acceptor_residue;
    bool acceptor_backbone = false;
    size_t acceptor_free_rot_pre;
    size_t acceptor_free_rot_bound;
    PROT::Atom* donor_atom;
    PROT::Atom* acceptor_atom;
    float distance;
    // equal operator that compares the donor and acceptor residues and backbones
    bool operator==(const PROT::SaltBridge& other){
        return (donor_residue->name() == other.donor_residue->name() && 
                acceptor_residue->name() == other.acceptor_residue->name() &&
                donor_backbone == other.donor_backbone && 
                acceptor_backbone == other.acceptor_backbone);
    }
};

// End the header guard from the start of the file
#endif
