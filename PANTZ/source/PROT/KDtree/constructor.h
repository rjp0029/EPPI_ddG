/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included by the KDtree.h header file. It implements
 * the constructor of the KDtree class. */

// Confirm that the KDtree class is loading the content
#ifndef KDtree_Loading_Status
#error KDtree methods must be included by KDtree.h
#endif

// The constructor of the KDtree class that instantiates the class
template <typename T>
PROT::KDtree<T>::KDtree () {
    this->root = nullptr;
}

// constructor for kdtree class with vector of template type
template <typename T>
PROT::KDtree<T>::KDtree (vector<T*> items) {
    this->m_item_ptrs = items;
    this->root = build_tree(items, 0);
}

// destructor for the kdtree class (i know its called constructor.h leave me alone)
template <typename T>
PROT::KDtree<T>::~KDtree () {
    delete_subtree(this->root);
}
