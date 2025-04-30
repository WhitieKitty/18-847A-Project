// EigenSolver.h
#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include "base_matrix.h"
#include <vector>

class EigenSolver {
public:
    // Compute the dominant eigenvalue and eigenvector using Power Iteration
    static std::pair<double, std::vector<double>> 
    powerIteration(const BaseMatrix& A, int maxIters = 1000, double tol = 1e-6);

    // Compute an eigenvalue near a target using Inverse Iteration
    static std::pair<double, std::vector<double>> 
    inverseIteration(const BaseMatrix& A, int maxIters = 1000, double tol = 1e-6);

    // Lanczos Algorithm for symmetric matrices
    static std::vector<double> lanczos(const BaseMatrix& A, int k, BaseMatrix& T, BaseMatrix& Q,int max_iters=100);

    // Arnoldi Iteration for general matrices
    static std::vector<double> arnoldi(const BaseMatrix& A, int k, BaseMatrix& H, BaseMatrix& Q,int max_iters=100);
    
    // QR Iteration for eigenvalues
    static std::vector<double> QRIteration(const BaseMatrix& A, int max_iters = 100, double tol = 1e-10);

};

#endif 
