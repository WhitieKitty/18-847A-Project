#include "Decomposition.h"
#include <cassert>
#include <cmath>

void Decomposition::LU(const BaseMatrix& A, BaseMatrix& L, BaseMatrix& U) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());
    for (int i = 0; i < n; ++i) {
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                double lij = L.get(i, j);
                double ujk = U.get(j, k);
                if (std::abs(lij) > 1e-12 && std::abs(ujk) > 1e-12)
                    sum += lij * ujk;
            }
            double aik = A.get(i, k);
            if (std::abs(aik) > 1e-12 || std::abs(sum) > 1e-12)
                U.set(i, k, aik - sum);
        }
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L.set(i, i, 1.0);
            } else {
                double sum = 0.0;
                for (int j = 0; j < i; ++j) {
                    double lkj = L.get(k, j);
                    double uji = U.get(j, i);
                    if (std::abs(lkj) > 1e-12 && std::abs(uji) > 1e-12)
                        sum += lkj * uji;
                }
                double aki = A.get(k, i);
                if (std::abs(aki) > 1e-12 || std::abs(sum) > 1e-12)
                    L.set(k, i, (aki - sum) / U.get(i, i));
            }
        }
    }
}

std::vector<double> Decomposition::solveLU(const BaseMatrix& L, const BaseMatrix& U, const std::vector<double>& b) {
    int n = L.getRows();
    assert(L.getRows() == L.getCols());
    assert(U.getRows() == U.getCols());
    assert((int)b.size() == n);
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            double lij = L.get(i, j);
            if (std::abs(lij) > 1e-12)
                y[i] -= lij * y[j];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            double uij = U.get(i, j);
            if (std::abs(uij) > 1e-12)
                x[i] -= uij * x[j];
        }
        x[i] /= U.get(i, i);
    }
    return x;
}

void Decomposition::QR(const BaseMatrix& A, BaseMatrix& Q, BaseMatrix& R) {
    int n = A.getRows();
    assert(A.getCols() == n);
    std::vector<std::vector<double>> q_vectors(n, std::vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i)
            q_vectors[k][i] = A.get(i, k);
        for (int j = 0; j < k; ++j) {
            double proj = 0.0;
            for (int l = 0; l < n; ++l) {
                double qjl = q_vectors[j][l];
                double alk = A.get(l, k);
                if (std::abs(qjl) > 1e-12 && std::abs(alk) > 1e-12)
                    proj += qjl * alk;
            }
            R.set(j, k, proj);
            for (int l = 0; l < n; ++l)
                q_vectors[k][l] -= proj * q_vectors[j][l];
        }
        double norm = 0.0;
        for (int i = 0; i < n; ++i)
            norm += q_vectors[k][i] * q_vectors[k][i];
        norm = std::sqrt(norm);
        R.set(k, k, norm);
        for (int i = 0; i < n; ++i)
            q_vectors[k][i] /= norm;
    }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            Q.set(i, j, q_vectors[j][i]);
}

void Decomposition::Cholesky(const BaseMatrix& A, BaseMatrix& L) {
    int n = A.getRows();
    assert(A.getCols() == n); // A must be square

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += L.get(i, k) * L.get(j, k);
            }
            if (i == j) {
                double value = A.get(i, i) - sum;
                assert(value > 0); // A must be positive definite
                L.set(i, j, std::sqrt(value));
            } else {
                L.set(i, j, (A.get(i, j) - sum) / L.get(j, j));
            }
        }
    }
}
