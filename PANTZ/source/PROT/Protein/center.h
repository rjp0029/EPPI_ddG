/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the Protein.h
 * header file, and includes pre-processor directives to that effect. It defines
 * the center method of the Protein class. */

// Make sure that the Protein class is currently being loaded
#ifndef ProteinClass_Loading_Status
#error Protein methods must be included by the Protein.h header file
#endif

// Move a Protein so that its center of mass is at the origin
PROT::Matrix PROT::Protein::center (const string how = "all") {
    // If the Protein is empty, throw an error
    if (m_count == 0) {
        string error = "An empty Protein cannot be centered.\n";
        throw PANTZ_error (error);}
    // Start by selecting the designated atoms
    vector<Atom *> ptrs; select_atoms(ptrs, how);
    // If no atoms were identified, throw an error
    if (ptrs.size() == 0) {
        string error = "The '" + how + "' atom selection criteria failed to "
                       "find any Atoms for Protein centering.\n";
        throw PANTZ_error (error);}
    // Create a 1 x Atom Coordinates matrix of zeros
    Matrix matrix (1, AtomCoordinates);
    // For each coordinate position
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Store the average value here
        coor average = 0.0;
        // Loop through the Atoms
        for(size_t j=0; j<ptrs.size(); ++j) {
            average += ptrs[j]->operator[](i);}
        matrix.set(0, i, average / ((coor) ptrs.size()));}
    // Move the Protein using those average coordinates
    private_move(&matrix, '-');
    return matrix;
}
