// Decomposition.h
#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include "base_matrix.h"
#include <vector>

class Decomposition {
public:
    static void LU(const BaseMatrix& A, BaseMatrix& L, BaseMatrix& U);

    static void QR(const BaseMatrix& A, BaseMatrix& Q, BaseMatrix& R);

    static void Cholesky(const BaseMatrix& A, BaseMatrix& L);

    static std::vector<double> solveLU(const BaseMatrix& L, const BaseMatrix& U, const std::vector<double>& b);
};

#endif 
