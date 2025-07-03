/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the CHARMM histidine correction methods of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Change HIS to HSD
void PROT::Protein::for_charmm_histidine_fix () {
    // If the protein is empty, be done
    if (m_count == 0) {return;}
    // Go through the residues
    for (size_t i=0; i<m_count; ++i) {
        // Use the corresponding method of the Residue class
        m_residues[i].for_charmm_histidine_fix();}
}

// Change HSD and HSE to HIS
void PROT::Protein::from_charmm_histidine_fix () {
    if (m_count == 0) {return;}
    for(size_t i=0; i<m_count; ++i) {
        m_residues[i].from_charmm_histidine_fix();}
}
