/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file is one that SHOULD BE EDITED by everyone who installs the PANTZ
 * code on their own system. It contains paths to key locations in your file
 * system, and that is information that is unique to each user. */

// A header guard to prevent this file from being included in a compiled program
// multiple times
#ifndef PANTZ_Macro_Guard
#define PANTZ_Macro_Guard 1

// The paths are stored as preprocessor macros. You should edit the information
// inside of the quotation marks.

// This is where you have installed the PANTZ code. This is likely a folder
// named PANTZ and it should contain the PANTZ.cpp program as well as the
// compiled executable after it has been compiled.
#define PANTZ_PATH "/Users/acr0116/Desktop/research/EPPI_ddg_project/EPPI_ddg/PANTZ/"

// If you want to use Rosetta through the PANTZ code, you must separately
// install it on your system. PANTZ works with three Rosetta programs: the
// minimization program (MIN), the interface analyzer (RIA), and the residue 
// energy breakdown protocol (REB). These macros should be the locations of the 
// corresponding executables.
//
// If you are not a member of the PROTEIN PANT(z) lab and you are accessing this
// file, it is almost certain that you need to use Rosetta to achieve the
// functionality you are looking for and therefore have to edit these paths.
#define ROSETTA_MIN_exec "/Users/acr0116/Desktop/rosetta/source/bin/minimize.default.macosclangrelease"
#define ROSETTA_RIA_exec "/Users/acr0116/Desktop/rosetta/source/bin/InterfaceAnalyzer.default.macosclangrelease"
#define ROSETTA_REB_exec "/Users/acr0116/Desktop/rosetta/source/bin/residue_energy_breakdown.default.macosclangrelease"

// Although they are not part of our publicly disseminated programs, the PROTEIN
// PANT(z) lab also does calculations with RosettaFold and CHARMM. Therefore, we
// have macros to those programs, too. If you are not a member of the PROTEIN
// PANT(z) lab, it is unlikely that you will use code that requires these macros
// and therefore you likely do not need to edit them.
#define ROSETTAFOLD_PATH "/path/to/rosettafold_repo"
#define CHARMM_exec "/path/to/charmm_executable"

// End the header guard from the start of the file
#endif
