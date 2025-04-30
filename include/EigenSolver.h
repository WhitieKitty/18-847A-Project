// EigenSolver.h
#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include "base_matrix.h"
#include <vector>

/**
 * @class EigenSolver
 * @brief Provides algorithms to compute eigenvalues and eigenvectors of matrices.
 */
class EigenSolver {
public:
    /**
     * @brief Computes the dominant eigenvalue and its corresponding eigenvector using the Power Iteration method.
     * 
     * @param A The input square matrix.
     * @param maxIters Maximum number of iterations (default: 1000).
     * @param tol Convergence tolerance (default: 1e-6).
     * 
     * @return A pair containing the dominant eigenvalue and the corresponding eigenvector.
     * 
     * @note A must be a square matrix.
     */
    static std::pair<double, std::vector<double>> 
    powerIteration(const BaseMatrix& A, int maxIters = 1000, double tol = 1e-6);

    /**
     * @brief Computes an eigenvalue and eigenvector near a given target using the Inverse Iteration method.
     * 
     * @param A The input square matrix.
     * @param maxIters Maximum number of iterations (default: 1000).
     * @param tol Convergence tolerance (default: 1e-6).
     * 
     * @return A pair containing the approximated eigenvalue and its corresponding eigenvector.
     * 
     * @note A must be invertible and square. The algorithm is sensitive to initial conditions.
     */
    static std::pair<double, std::vector<double>> 
    inverseIteration(const BaseMatrix& A, int maxIters = 1000, double tol = 1e-6);

    /**
     * @brief Performs the Lanczos algorithm to approximate eigenvalues of a symmetric matrix.
     * 
     * @param A The input symmetric matrix.
     * @param k Number of steps (dimension of Krylov subspace).
     * @param T The resulting tridiagonal matrix.
     * @param Q The orthonormal basis of the Krylov subspace.
     * @param max_iters Maximum number of iterations (default: 100).
     * 
     * @return A vector of approximate eigenvalues (diagonal of T).
     * 
     * @note A must be symmetric and square.
     */
    static std::vector<double> lanczos(const BaseMatrix& A, int k, BaseMatrix& T, BaseMatrix& Q, int max_iters = 100);

    /**
     * @brief Performs the Arnoldi iteration to approximate eigenvalues of a general (possibly non-symmetric) matrix.
     * 
     * @param A The input square matrix.
     * @param k Number of steps (dimension of Krylov subspace).
     * @param H The resulting upper Hessenberg matrix.
     * @param Q The orthonormal basis of the Krylov subspace.
     * @param max_iters Maximum number of iterations (default: 100).
     * 
     * @return A vector of approximate eigenvalues (eigenvalues of H).
     */
    static std::vector<double> arnoldi(const BaseMatrix& A, int k, BaseMatrix& H, BaseMatrix& Q, int max_iters = 100);

    /**
     * @brief Computes all eigenvalues of a matrix using the QR Iteration method.
     * 
     * @param A The input square matrix.
     * @param max_iters Maximum number of iterations (default: 100).
     * @param tol Convergence tolerance (default: 1e-10).
     * 
     * @return A vector of approximated eigenvalues.
     * 
     * @note The method does not return eigenvectors. A must be square.
     */
    static std::vector<double> QRIteration(const BaseMatrix& A, int max_iters = 100, double tol = 1e-10);
};

#endif // EIGENSOLVER_H
