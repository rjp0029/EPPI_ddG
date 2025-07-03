/* Created by the clay at Auburn University.
 *
 * This file contains the declaration of the KDtree class for PDB calculations.
 * It also includes the header files where the methods of the class are
 * implemented. */

// Use a header guard to make sure this file is only included in a compiled
// program a single time
#ifndef Proteins_KDtree_Guard
#define Proteins_KDtree_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error KDtree.h must be included by Proteins.h
#endif

// Define the KDtree class
template <class T>
class PROT::KDtree {
    private:
        // the node struct which contains the data and the left and right children
        struct Node {
            T* item;
            Node* left;
            Node* right;

            Node(T* d) : item(d), left(nullptr), right(nullptr) {}
        };
        size_t dimension; // Dimension to split on: 0 for x, 1 for y, 2 for z
        // function to build the tree
        KDtree::Node * build_tree(vector<T*> items, size_t depth);
        // item pointers for the tree
        vector<T*> m_item_ptrs;
        // the root node
        Node* root;
        // the radius search function
        void radius_search(Node* node, T* data, float radius, vector<T*>& neighbors, size_t depth=0);
    public:
        // default constructor
        KDtree();
        // constructor from a vector of items
        KDtree(vector<T*> data);
        // destructor
        ~KDtree();
        // delete subtree
        void delete_subtree(Node* node);
        // find the nearest neighbor
        vector<T*> nearest_neighbors(T* data, size_t num_neighbors);
        // radius neighbors
        vector<T*> radius_neighbors(T* data, float radius);
        // get the items in the tree
        void get_items(Node* node, vector<T*>& data);

    // End the class definition
};

// Define a preprocessor variable to guarantee that the KDtree methods are
// included here and only here
#define KDtree_Loading_Status 1

// Include the files that implement class methods
#include "KDtree/constructor.h"
#include "KDtree/build_tree.h"
#include "KDtree/nearest_neighbors.h"
#include "KDtree/radius_neighbors.h"

// Undefine the matrix loading status variable
#undef KDtree_Loading_Status

// End the header guard from the start of the file
#endif
