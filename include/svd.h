#ifndef SVD_H
#define SVD_H

#include "base_matrix.h"
#include <vector>

/**
 * @class SVD
 * @brief Provides algorithms for computing the largest singular value decomposition (SVD) of a matrix.
 */
class SVD {
public:
    /**
     * @brief Computes the largest singular value of matrix A using the Lanczos Bidiagonalization algorithm.
     * 
     * @param A The input matrix.
     * @param k Number of Lanczos steps (i.e., the size of the bidiagonal matrix to compute).
     * 
     * @return The largest singular value approximated from the bidiagonal matrix.
     * 
     * @note A can be any real matrix (not necessarily square). The accuracy depends on k.
     */
    static double lanczos_bidiag(const BaseMatrix& A, int k);
};

#endif // SVD_H
