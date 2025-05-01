#define ACCELERATE_NEW_LAPACK
#include "svd.h"
#include <cassert>
#include <cmath>
#include <random>
#include <vector>
#include <numeric>
#include <iostream>
#include <Accelerate/Accelerate.h>

double SVD::lanczos_bidiag(const BaseMatrix& A, int k) {
    std::vector<std::vector<double>> V(k+1, std::vector<double>(A.getCols()));
    std::vector<std::vector<double>> U(k,   std::vector<double>(A.getRows()));
    std::vector<double> alpha(k), beta(k+1);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < V[0].size(); i++) {
        V[0][i] = dis(gen);
    }

    double norm = std::inner_product(V[0].begin(), V[0].end(), V[0].begin(), 0.0);
    norm = std::sqrt(norm);
    for (int i = 0; i < V[0].size(); i++) {
        V[0][i] /= norm;
    }

    const BaseMatrix* AT = A.transpose();

    for (int i = 0; i < k; i++) {
        std::vector<double> w = A.multiply(V[i]);
        if (i > 0) {
            for (int j = 0; j < (int)w.size(); j++)
                w[j] -= beta[i] * U[i-1][j];
        }

        alpha[i] = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
        for (int j = 0; j < (int)w.size(); j++)
            U[i][j] = w[j] / alpha[i];

        std::vector<double> z = AT->multiply(U[i]);
        for (int j = 0; j < (int)z.size(); j++)
            z[j] -= alpha[i] * V[i][j];

        beta[i+1] = std::sqrt(std::inner_product(z.begin(), z.end(), z.begin(), 0.0));
        for (int j = 0; j < (int)z.size(); j++)
            V[i+1][j] = z[j] / beta[i+1];
    }

    int m = k, n = k+1, lda = m;
    int ldu = m, ldvt = n;
    std::vector<double> B(m * n, 0.0);
    // construct Bidiagonal Matrix
    for (int i = 0; i < k; ++i) {
        B[i*n + i] = alpha[i];
        B[i*n + i + 1] = beta[i+1];
    }

    int info;
    int mn = std::min(m,n), mx = std::max(m,n);
    int lwork = -1;
    double wkopt;
    double* work;
    double s[n], u[ldu*m], vt[ldvt*n];
    // get optimal workspace size
    dgesvd_( "All", "All", &m, &n, B.data(), &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    // compute SVD
    dgesvd_(
    "All", "All", &m, &n,
    B.data(), &lda,
    s, u, &ldu,
    vt, &ldvt,
    work, &lwork,
    &info
    );
    if (info != 0) {
        std::cerr << "dgesvd_ failed, info=" << info << "\n";
        delete AT;
        return -1;
    }

    delete AT;
    std::cout << std::endl;
    return s[0];
}