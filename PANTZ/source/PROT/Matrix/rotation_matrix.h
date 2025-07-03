/* Created by the Pantazes Lab at Auburn University.
 *
 * This file is intended to be included in a compiled program by the main
 * Matrix.h header file for this class. The file contains the code to calculate the 
 * rotation matrix that superimposes this matrix (set of coordinates) onto the reference
 * (i.e., part of the kabsch algorithm to align proteins with minimal RMSD). */

// Error check the inclusion chain for the file
#ifndef Matrix_Loading_Status
#error Matrix methods must be included by the Matrix.h header file
#endif

#include "../../external/Eigen/Dense"

// singular value decomposition to to get the rotation matrix to put one on another
PROT::Matrix PROT::Matrix::rotation_matrix (const Matrix * other) const {
    // if the matrices are not the same dimensions, throw an error
    if (m_rows != other->m_rows || m_columns != other->m_columns) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        stringstream c3; c3 << other->m_rows;
        stringstream c4; c4 << other->m_columns;
        string error = "It is not possible to calculate the rotation matrix of a "
                     + c1.str() + "x" + c2.str() + " matrix and a "
                     + c3.str() + "x" + c4.str() + " matrix.\n";
        throw PANTZ_error (error);}
    // turn this matrix into an Eigen matrix
    Eigen::MatrixXf target(m_rows, m_columns);
    for (size_t i=0; i<m_rows; ++i) {
        for (size_t j=0; j<m_columns; ++j) {
            target(i, j) = operator()(i, j);
        }
    }
    // turn the other matrix into an Eigen matrix
    Eigen::MatrixXf reference(other->m_rows, other->m_columns);
    for (size_t i=0; i<other->m_rows; ++i) {
        for (size_t j=0; j<other->m_columns; ++j) {
            reference(i, j) = other->operator()(i, j);
        }
    }
    // get the centroid of the reference protein
    Eigen::Vector3f ref_centroid = reference.colwise().mean();
    // get the centroid of the target protein
    Eigen::Vector3f target_centroid = target.colwise().mean();
    // translate the proteins to the origin (for calculation only, doesnt change the protein)
    reference.rowwise() -= ref_centroid.transpose();
    target.rowwise() -= target_centroid.transpose();
    // calculate the covariance matrix
    Eigen::MatrixXf covariance = target.transpose() * reference;
    // calculate the singular value decomposition
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(covariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // calculate the rotation matrix
    Eigen::MatrixXf rotation = svd.matrixU() * svd.matrixV().transpose();
    // turn the rotation matrix into a pantz matrix
    PROT::Matrix rot_matrix(rotation.rows(), rotation.cols());
    for (size_t i=0; i<rotation.rows(); ++i) {
        for (size_t j=0; j<rotation.cols(); ++j) {
            rot_matrix.set(i, j, rotation(i, j));
        }
    }
    return rot_matrix;
}
