/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the declaration of the NeighborBox class, which is used to
 * rapidly find nearby atoms for various calculations. This file is responsible
 * for including the header files of the class */

// Use a header guard to make sure the file is only included in a compiled
// program a single time
#ifndef Proteins_Neighbor_Guard
#define Proteins_Neighbor_Guard 1

// Make sure this file is being included in the manner that is expected
#ifndef Proteins_Loading_Status
#error NeighborBox.h must be included by Proteins.h
#endif

// Define the class
class PROT::NeighborBox {

    // The information stored in the class is private
    private:
        // The minimum and maximum X, Y, and Z coordinates
        coor m_XL, m_XU, m_YL, m_YU, m_ZL, m_ZU;
        // The length of a bin side
        coor m_bin_side;
        // The number of bins in the X, Y and Z directions
        size_t m_Xbins, m_Ybins, m_Zbins;
        // The stored pointers to the atoms that make up the system
        vector<vector<vector<vector<Atom *> > > > m_atom_pointers;
        // The number of stored atoms
        size_t m_size;

    // Private methods to control class behavior
    private:
        // A method to assign default values to the class variables
        void initialize ();
        // Copy content from one instance of the class to another
        void copy (const NeighborBox&);
        // Calculate the bin of a point
        size_t calculate_bin (const coor, const char);
        // Construct a NeighborBox object
        void construct (vector<Atom *>&, const coor, const bool);

    // Public methods of the class
    public:
        // A default class constructor
        NeighborBox () {initialize();}
        // A series of standard class constructors
        NeighborBox (vector<Atom>&, const coor, const bool);
        NeighborBox (vector<AtomPtr>&, const coor, const bool);
        NeighborBox (vector<Residue>&, const coor, const bool);
        NeighborBox (vector<ResiduePtr>&, const coor, const bool);
        NeighborBox (vector<Protein>&, const coor, const bool);
        NeighborBox (vector<Protein>&, vector<Protein>&, const coor, const bool);
        // Copy construction and assignment
        NeighborBox (const NeighborBox& other) {initialize(); copy(other);}
        void operator= (const NeighborBox& other) {copy(other);}
        // The number of atoms
        size_t size () const {return m_size;}
        // Find candidate neighbors to a point
        void find_neighbors (vector<Atom *>&, const coor, const coor, 
                             const coor, const coor);
        void find_neighbors (vector<AtomPtr>&, const coor, const coor,
                             const coor, const coor);

    // End the class definition
};

// Include the header files of the class
#define NeighborBox_Loading_Status 1
#include "Neighbor/initialize.h"
#include "Neighbor/copy.h"
#include "Neighbor/calculate.h"
#include "Neighbor/construct.h"
#include "Neighbor/find.h"

// End the loading status preprocessor directive
#undef NeighborBox_Loading_Status

// End the header guard from the start of the file
#endif
