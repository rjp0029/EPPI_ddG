/* Created by the PROTEIN PANT(z) Lab at Auburn University. 
 *
 * This file declares a namespace for our "expected persistent pairwise
 * interation" (EPPI) features, which were initially published in doi:
 * 10.1002/prot.26773 ("Using Short Molecular Dynamics Simulations to Determine
 * the Important Features of Interactions in Antibody-Protein Complexes"). */

// Use a header guard to make sure the file is only included in a compiled
// program a single time
#ifndef EPPI_header_guard
#define EPPI_header_guard 1

// This file should be included in PANTZ calculations by Methods.h
#ifndef Methods_Loading_Status
#error EPPI.h must be included by Methods.h
#endif

// Declare the namespace
namespace EPPI {

    // The names of the base EPPI features (i.e., the ones that the code
    // actually directly calculates) are:
    const string BaseFeatureNames = "n_epp_sb" // 0: number (n) of expected persistent pairwise (epp) salt bridges (sb)
        " n_epp_hb"     // 1: n of epp hydrogen bonds (hb)
        " n_epp_hp"     // 2: n of epp hydrophbic interactions (hp)
        " n_np_sb"      // 3: n of non-persistent (np) sb
        " n_np_hb"      // 4: n of np hb
        " n_np_hp"      // 5: n of np hp
        " rl_epp_sb"    // 6: rotamers lost (rl) of epp sb
        " rl_epp_hb"    // 7: rl of epp hb
        " rl_epp_hp"    // 8: rl of epp hp
        " rl_np_sb"     // 9: rl of np sb
        " rl_np_hb"     // 10: rl of np hb
        " rl_np_hp"     // 11: rl of np hp
        " dg_separated" // 12: Rosetta-calculated dG_separated
        " bsa"          // 13: Buried surface area as calculated by Rosetta
        " sc"           // 14: Rosetta-calculated shape complementarity
        " epbe_sb"      // 15: Expected persistent binding energy (epbe) of sb
        " epbe_hb"      // 16: epbe of hb
        " epbe_hp"      // 17: epbe of hp
        " n_st_st_sb"   // 18: n of stable (st) - st sb
        " n_st_st_hb"   // 19: n of st - st hb
        " n_st_st_hp";  // 20: n of st - st hp

    // Split them into a vector
    const vector<string> BaseFeaturesList = Text::split(BaseFeatureNames);

    // A string of the primary features used in the EPPI calculations. These are
    // the base features, as well as certain linear combinations of them.
    const string PrimaryFeatureNames = BaseFeatureNames +
        " n_epp_polar"  // 21: n_epp_sb + n_epp_hb
        " n_eppi"       // 22: n of epp interactions (eppi)  = n_epp_polar + n_epp_hp
        " n_np_polar"   // 23: n_np_sb + n_np_hb
        " n_nppi"       // 24: n_np_polar + n_np_hp
        " rl_epp_polar" // 25: rl_epp_sb + rl_epp_hb
        " rl_eppi"      // 26: rl_epp_polar + rl_epp_hp
        " rl_np_polar"  // 27: rl_np_sb + rl_np_hb
        " rl_nppi"      // 28: rl_np_polar + rl_np_hp
        " epbe_polar"   // 29: epbe_sb + epbe_hb
        " epbe"         // 30: epbe_polar + epbe_hp
        " npbe"         // 31: non-persistent binding energy = dg_separated - epbe
        " n_st_st_polar"// 32: n_st_st_sb + n_st_st_hb
        " n_st_st_pi"   // 33: n_st_st_polar + n_st_st_hp
        " n_total"      // 34: n_eppi + n_nppi
        " rl_total";    // 35: rl_eppi + rl_nppi

    // Split them into a vector
    const vector<string> PrimaryFeaturesList = Text::split(PrimaryFeatureNames);

    // A function to scan the interface of a binding complex of proteins and get
    // the expected persistent pairwise interaction features. It writes them to
    // a file, features.txt, in the output directory and also runs the Rosetta
    // interface analyzer and per residue analysis
    void calculate_eppi_features (vector<PROT::Protein *>&, string&, string&, bool);

    // Gather the EPPI features from directories containing the required files
    void gather_eppi_features (string, vector<string>, string, string, bool);

    // Classes that are specific to the EPPI calculations. They are distinct
    // from those in the PROT namespace and are much lighter.
    class Residue;
    class Interaction;
    class HydrogenBond;
    class SaltBridge;
    class Hydrophobic;

    // The Feature class is primarily naming information about features,
    // including their name, index in a list of features being used in an
    // analysis, and how it is calculated from the primary features
    class Feature;

    // The Binding Complex Properties (BCProps) class are the actual calculated
    // feature values for particular binding complex. It starts from the base
    // features, calculates the primary features, and then uses the primary
    // features to populate a vector of values for the specific features that
    // are used in an analysis.
    class BCProps;

    // End the namespace
};

// Include EPPI namespace files

// Define a loading status flag
#define EPPI_Loading_Status 1

// Include the header files
#include "EPPI/Feature.h"
#include "EPPI/calculate_eppi_features.h"
#include "EPPI/Residue.h"
#include "EPPI/Interaction.h"
#include "EPPI/HydrogenBond.h"
#include "EPPI/SaltBridge.h"
#include "EPPI/Hydrophobic.h"
#include "EPPI/BCProps.h"
#include "EPPI/gather_eppi_features.h"

// Undefine the loading status
#undef EPPI_Loading_Status

// End the header guard from the start of the file
#endif
