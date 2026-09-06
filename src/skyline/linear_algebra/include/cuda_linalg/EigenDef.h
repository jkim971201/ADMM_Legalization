#ifndef EIGEN_DEF_H
#define EIGEN_DEF_H

// Eigen
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

namespace cuda_linalg 
{

typedef Eigen::Triplet<double>                       EigenTriplet;
typedef Eigen::VectorXd                              EigenVector;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSMatrix; // CSR (Default : CSC)
typedef Eigen::MatrixXd                              EigenDMatrix;
typedef Eigen::SelfAdjointEigenSolver<EigenDMatrix>  EigenSolver;

}

#endif
