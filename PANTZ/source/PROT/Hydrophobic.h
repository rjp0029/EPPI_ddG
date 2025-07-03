/* Created by the clay at Auburn University.
 *
 * This file contains the declaration of the Hydrophobic struct. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_Hydrophobic_Guard
#define Proteins_Hydrophobic_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error Hydrophobic.h must be included by Proteins.h
#endif

// confirm that the Residue header files have been included
#ifndef Proteins_Residue_Guard
#error Hydrophobic.h must be included after Residue.h
#endif

struct PROT::Hydrophobic{
    PROT::Residue* res1;
    bool res1_backbone = false;
    size_t res1_free_rot_pre;
    size_t res1_free_rot_bound;
    PROT::Residue* res2;
    bool res2_backbone = false;
    size_t res2_free_rot_pre;
    size_t res2_free_rot_bound;
    float phobic_BSASA;
    // equal operator that compares the donor and acceptor residues and backbones
    bool operator==(const PROT::Hydrophobic& other){
        return (res1 == other.res1 && 
                res2 == other.res2 &&
                res1_backbone == other.res1_backbone && 
                res2_backbone == other.res2_backbone);
    }
};

// End the header guard from the start of the file
#endif
