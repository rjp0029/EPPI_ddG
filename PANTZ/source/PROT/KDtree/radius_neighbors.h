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
vector<T*> PROT::KDtree<T>::radius_neighbors(T* item_ptr, float radius) {
    // handle null pointer
    if (root == nullptr) {
        return vector<T*>();
    }
    // get the neighbors
    vector<T*> neighbors;
    radius_search(root, item_ptr, radius, neighbors);
    // sort the neighbors by distance to the item
    sort(neighbors.begin(), neighbors.end(), [item_ptr](T* a, T* b) {
        return a->distance(*item_ptr) < b->distance(*item_ptr);
    });
    // remove this item from the neighbors
    neighbors.erase(remove(neighbors.begin(), neighbors.end(), item_ptr), neighbors.end());
    // return the neighbors
    return neighbors;
}

// // radius search
// template <typename T>
// void PROT::KDtree<T>::radius_search(Node* node, T* data, float radius, vector<T*>& neighbors, size_t depth) {
//     // handle null pointer
//     if (node == nullptr) {
//         return;
//     }
//     // get the dimension
//     size_t dim = depth % 3;
//     // get the distance between the node and the data
//     float dist = node->item->distance(*data);
//     // if the distance is less than the radius, add the node to the neighbors
//     if (dist < radius) {
//         neighbors.push_back(node->item);
//     }
//     // if dim is 0, check the x coordinate, if 1, check the y coordinate, if 2, check the z coordinate
//     if (dim == 0) {
//         // if the distance is less than the radius, search the left and right children
//         if (node->left != nullptr and node->item->x() - data->x() > -radius) {
//             radius_search(node->left, data, radius, neighbors, depth + 1);
//         }
//         if (node->right != nullptr and node->item->x() - data->x() < radius) {
//             radius_search(node->right, data, radius, neighbors, depth + 1);
//         }
//     } else if (dim == 1) {
//         // if the distance is less than the radius, search the left and right children
//         if (node->left != nullptr and node->item->y() - data->y() > -radius) {
//             radius_search(node->left, data, radius, neighbors, depth + 1);
//         }
//         if (node->right != nullptr and node->item->y() - data->y() < radius) {
//             radius_search(node->right, data, radius, neighbors, depth + 1);
//         }
//     } else {
//         // if the distance is less than the radius, search the left and right children
//         if (node->left != nullptr and node->item->z() - data->z() > -radius) {
//             radius_search(node->left, data, radius, neighbors, depth + 1);
//         }
//         if (node->right != nullptr and node->item->z() - data->z() < radius) {
//             radius_search(node->right, data, radius, neighbors, depth + 1);
//         }
//     }
// }

template <typename T>
void PROT::KDtree<T>::radius_search(Node* node, T* data, float radius, vector<T*>& neighbors, size_t depth) {
    // handle null pointer
    if (node == nullptr) {
        return;
    }

    // get the dimension
    size_t dim = depth % 3;

    // get the distance between the node and the data
    float dist = node->item->distance(*data);

    // if the distance is less than the radius, add the node to the neighbors
    if (dist < radius) {
        neighbors.push_back(node->item);
    }

    // check which side of the splitting plane the data point lies on
    float diff;
    if (dim == 0) {
        diff = data->x() - node->item->x();
    } else if (dim == 1) {
        diff = data->y() - node->item->y();
    } else {
        diff = data->z() - node->item->z();
    }

    if (diff <= 0) {
        // search the left subtree
        radius_search(node->left, data, radius, neighbors, depth + 1);
        // if there could be points on the right within the radius, search the right subtree
        if (diff * diff < radius * radius) {
            radius_search(node->right, data, radius, neighbors, depth + 1);
        }
    } else {
        // search the right subtree
        radius_search(node->right, data, radius, neighbors, depth + 1);
        // if there could be points on the left within the radius, search the left subtree
        if (diff * diff < radius * radius) {
            radius_search(node->left, data, radius, neighbors, depth + 1);
        }
    }
}
