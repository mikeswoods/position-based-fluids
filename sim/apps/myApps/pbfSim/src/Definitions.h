/*******************************************************************************
 * Definitions.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// For the Eigen library:
#include "Core"
#include "Dense"
#include "Sparse"

// Basic Eigen vectors and matrices shorthand types:
// This was adapted from code by Tiantian Liu in the deformable mesh assignment

typedef Eigen::Matrix<float, 3, 3, 0, 3 ,3> EigenMatrix3;
typedef Eigen::Matrix<float, 3, 1, 0, 3 ,1> EigenVector3;
typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::Triplet<float, int> SparseMatrixTriplet;

/******************************************************************************/

#define block_vector(a) block<3,1>(3*(a), 0)

/******************************************************************************/

#endif