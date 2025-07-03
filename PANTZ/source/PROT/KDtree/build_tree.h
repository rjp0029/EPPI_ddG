/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the KDtree.h header file. It 
    * implements the function to build the treee */

// Confirm that the KDtree class is loading the content
#ifndef KDtree_Loading_Status
#error KDtree methods must be included by KDtree.h
#endif

// build the tree for the KDtree class
template <typename T>
typename PROT::KDtree<T>::Node* PROT::KDtree<T>::build_tree(vector<T*> items, size_t depth) {
    if (items.empty()) {
        return nullptr;
    }

    size_t axis = depth % 3;

    // Sort items based on the dimension
    if (axis == 0) {
        sort(items.begin(), items.end(), [](T* a, T* b) {
            return a->x() < b->x();
        });
    } else if (axis == 1) {
        sort(items.begin(), items.end(), [](T* a, T* b) {
            return a->y() < b->y();
        });
    } else {
        sort(items.begin(), items.end(), [](T* a, T* b) {
            return a->z() < b->z();
        });
    }

    size_t median = items.size() / 2;
    PROT::KDtree<T>::Node* node = new PROT::KDtree<T>::Node(items[median]);
    // left and right children
    vector<T*> left_items(items.begin(), items.begin() + median);
    vector<T*> right_items(items.begin() + median + 1, items.end());
    node->left = build_tree(left_items, depth + 1);
    node->right = build_tree(right_items, depth + 1);
    return node;
}

// delete subtree for KDtree class with template type
template <typename T>
void PROT::KDtree<T>::delete_subtree(Node* node) {
    if (node == nullptr) {
        return;
    }
    delete_subtree(node->left);
    delete_subtree(node->right);
    delete node;
}
