// Decomposition.h
#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include "base_matrix.h"
#include <vector>

/**
 * @class Decomposition
 * @brief A utility class that provides matrix decomposition algorithms: LU, QR, and Cholesky.
 */
class Decomposition {
public:
    /**
     * @brief Performs LU decomposition on matrix A such that A = L * U.
     * 
     * @param A The input square matrix to decompose.
     * @param L The output lower triangular matrix.
     * @param U The output upper triangular matrix.
     * 
     * @note The function assumes A is square and invertible.
     */
    static void LU(const BaseMatrix& A, BaseMatrix& L, BaseMatrix& U);
    
    /**
     * @brief Performs QR decomposition using the Gram-Schmidt process.
     * 
     * @param A The input square matrix.
     * @param Q The output orthogonal matrix.
     * @param R The output upper triangular matrix.
     * 
     * @note The function assumes A is a square matrix.
     */
    static void QR(const BaseMatrix& A, BaseMatrix& Q, BaseMatrix& R);

    /**
     * @brief Performs Cholesky decomposition on matrix A such that A = L * L^T.
     * 
     * @param A The input symmetric positive definite matrix.
     * @param L The output lower triangular matrix.
     * 
     * @note The function assumes A is symmetric, positive definite, and square.
     */
    static void Cholesky(const BaseMatrix& A, BaseMatrix& L);

    /**
     * @brief Solves Ax = b using LU decomposition (A = L * U).
     * 
     * @param L The lower triangular matrix from LU decomposition.
     * @param U The upper triangular matrix from LU decomposition.
     * @param b The right-hand side vector.
     * 
     * @return A vector x that solves Ax = b.
     * 
     * @note The function uses forward and backward substitution.
     */
    static std::vector<double> solveLU(const BaseMatrix& L, const BaseMatrix& U, const std::vector<double>& b);
};

#endif // DECOMPOSITION_H
