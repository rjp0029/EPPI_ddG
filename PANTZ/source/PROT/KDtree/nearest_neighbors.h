/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the KDtree.h header file. It 
    * implements the function to find the k nearest neighbors to an item*/

// Confirm that the KDtree class is loading the content
#ifndef KDtree_Loading_Status
#error KDtree methods must be included by KDtree.h
#endif

// get the nearest neighbors
template <typename T>
vector<T*> PROT::KDtree<T>::nearest_neighbors(T* item_ptr, size_t num_neighbors) {
    // handle null pointer
    if (item_ptr == nullptr) {
        return vector<T*>();
    }
    vector<T*> neighbors;
    vector<T*> items;
    get_items(this->root, items);
    // sort the items based on the distance from the item_ptr
    sort(items.begin(), items.end(), [item_ptr](T* a, T* b) {
        return item_ptr->distance(*a) < item_ptr->distance(*b);
    });
    // get the first num_neighbors items
    for (size_t i = 0; i < num_neighbors; i++) {
        neighbors.push_back(items[i]);
    }
    return neighbors;
}

// get items for kdtree class with template type
template <typename T>
void PROT::KDtree<T>::get_items(Node* node, vector<T*>& items) {
    if (node == nullptr) {
        return;
    }
    items.push_back(node->item);
    get_items(node->left, items);
    get_items(node->right, items);
}