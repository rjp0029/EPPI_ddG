/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be loaded directly from the Residue.h header file,
 * and has preprocessor directives to control that behavior. It contains the
 * histidine naming correction methods of the Residue class. */

// Make sure that the Residue class is currently loading methods
#ifndef Residue_Loading_Status
#error Methods of the Residue class must be loaded from the Residue.h header file
#endif

// If the residue is going to CHARMM and is histidine, change it from HIS to HSD
void PROT::Residue::for_charmm_histidine_fix () {
    if (((m_name == "HIS") || (m_name == "HID")) || ((m_name == "HIE") || (m_name == "HIP"))) {
        m_name = "HSD";
        if (m_count > 0) {
            for(size_t i=0; i<m_count; ++i) {
                m_atoms[i].m_residue = "HSD";}}}
}

// If the Residue's name should be HIS, fix that.
void PROT::Residue::from_charmm_histidine_fix () {
    if ((m_name == "HSD") || (m_name == "HSE")) {
        m_name = "HIS";
        if (m_count > 0) {
            for(size_t i=0; i<m_count; ++i) {
                m_atoms[i].m_residue = "HIS";}}}
}
