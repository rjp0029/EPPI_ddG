/* Created by clay at Auburn University.
*
* This file implements the functions to set the m_sasa_points for 
* all atoms in a given vector of atoms. It makes a kd tree of the
* atoms and then sets the sasa points for each atom.
*
*/

// This file is supposed to be loaded by Methods.h
#ifndef Methods_Loading_Status
#error General Methods must be included by Methods.h
#endif

// include the KDtree class
#include "../PROT/KDtree.h"


// a function to set the sasa points for all atoms in a given vector of atoms 
// (the points on the atoms that are not buried by other atoms)
void METHODS::set_sasa_points(vector<PROT::Atom*> atoms){
    // create a KDtree of the atoms
    PROT::KDtree<PROT::Atom> kd_tree(atoms);
    // set the sasa points on the atoms in the 
    for (auto atom : atoms) {
        // get neighbors of the atom
        vector<PROT::Atom*> neighbors = kd_tree.radius_neighbors(atom, 8.0);
        // get the mesh points of the atom (with vdw radius)
        vector<vector<float>> mesh = atom->fibonacci_sphere(80, 1.4);
        // loop through the points in the mesh
        for (auto point : mesh) {
            // go through the atoms in the neighbors and see if the point is accessible
            bool solvent_accessible = true;
            // go through this atoms neighbors and determine if the point is buried by an atom 
            // on the same chain, if it is then it is not solvent accessible, otherwise it is
            for (auto atom2 : neighbors) {
                if ((atom == atom2) or (atom->protein() != atom2->protein())) {
                    continue;
                }
                // if the distance between the point and the atom is less than the sum of the radii
                float dist = sqrt(pow(atom2->x() - point[0], 2) + 
                                pow(atom2->y() - point[1], 2) + 
                                pow(atom2->z() - point[2], 2));
                if (dist < atom2->lj_sigma() + 1.4) {
                    solvent_accessible = false;
                    break;
                }
            }
            // if the point is solvent accessible, add it to the atom's sasa points
            if (solvent_accessible) {
                atom->add_sasa_point(point);
            }
        }
    }
}