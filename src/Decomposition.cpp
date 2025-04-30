#include "Decomposition.h"
#include <cassert>
#include <cmath>

// LU Decomposition
// This function decomposes a matrix A into L and U such that A = L * U
// L is a lower triangular matrix and U is an upper triangular matrix.
// The function assumes that A is a square matrix and that it is invertible.
void Decomposition::LU(const BaseMatrix& A, BaseMatrix& L, BaseMatrix& U) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols()); // A must be square

    for (int i = 0; i < n; ++i) {
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                // L[i][j] * U[j][k]
                // We need to check if lij and ujk are not too small to avoid numerical issues
                double lij = L.get(i, j);
                double ujk = U.get(j, k);
                if (std::abs(lij) > 1e-12 && std::abs(ujk) > 1e-12)
                    sum += lij * ujk;
            }
            double aik = A.get(i, k);
            // We need to check if aik and sum are not too small to avoid numerical issues
            if (std::abs(aik) > 1e-12 || std::abs(sum) > 1e-12)
                U.set(i, k, aik - sum);
        }
        for (int k = i; k < n; ++k) {
            // set diagonal elements of L to 1
            if (i == k) {
                L.set(i, i, 1.0);
            } else {
                // L[k][i] = (A[k][i] - sum) / U[i][i]
                double sum = 0.0;
                for (int j = 0; j < i; ++j) {
                    // Sum L[k][j] * U[j][i]
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

// This function solves the system of equations Ax = b using the LU decomposition
// where A = L * U. It first solves Ly = b for y using forward substitution
// and then solves Ux = y for x using backward substitution.
// The function assumes that L and U are the result of a LU decomposition of A
// and that b is a vector of the same size as A's rows.
// The function returns the solution vector x.
std::vector<double> Decomposition::solveLU(const BaseMatrix& L, const BaseMatrix& U, const std::vector<double>& b) {
    int n = L.getRows();
    assert(L.getRows() == L.getCols());
    assert(U.getRows() == U.getCols());
    assert((int)b.size() == n);
    std::vector<double> y(n, 0.0);
    // Forward substitution to solve Ly = b
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            double lij = L.get(i, j);
            // We need to check if lij is not too small to avoid numerical issues
            if (std::abs(lij) > 1e-12)
                y[i] -= lij * y[j];
        }
    }
    std::vector<double> x(n, 0.0);
    // Backward substitution to solve Ux = y
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

// QR Decomposition using Gram-Schmidt Process
// This function decomposes a square matrix A into an orthogonal matrix Q
// and an upper triangular matrix R such that A = Q * R.
void Decomposition::QR(const BaseMatrix& A, BaseMatrix& Q, BaseMatrix& R) {
    int n = A.getRows();
    assert(A.getCols() == n);
    // q_vectors will hold the orthogonal vectors (columns of Q)
    std::vector<std::vector<double>> q_vectors(n, std::vector<double>(n, 0.0));
    
    for (int k = 0; k < n; ++k) {
        // Initialize q_vectors[k] to the k-th column of A
        for (int i = 0; i < n; ++i) {
            q_vectors[k][i] = A.get(i, k);
        }
        // Subtract the projections of q_vectors[k] onto the previous q_vectors
        for (int j = 0; j < k; ++j) {
            double proj = 0.0;
            for (int l = 0; l < n; ++l) {
                double qjl = q_vectors[j][l];
                double alk = A.get(l, k);
                if (std::abs(qjl) > 1e-12 && std::abs(alk) > 1e-12) {
                    proj += qjl * alk;
                }
            }
            // R[j][k] = <q_j, a_k> (dot product)
            R.set(j, k, proj);
            for (int l = 0; l < n; ++l)
                q_vectors[k][l] -= proj * q_vectors[j][l];
        }
        // Normalize q_vectors[k] to get the k-th column of Q
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += q_vectors[k][i] * q_vectors[k][i];
        }
        norm = std::sqrt(norm);
        // R[k][k] = ||q_k|| (norm)
        R.set(k, k, norm);
        for (int i = 0; i < n; ++i) {
            q_vectors[k][i] /= norm;
        }
    }
    // Set the orthogonal matrix Q
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            Q.set(i, j, q_vectors[j][i]);
        }
    }
}

// Cholesky Decomposition
// This function decomposes a symmetric positive definite matrix A into
// a lower triangular matrix L such that A = L * L^T.
// The function assumes that A is a square matrix and that it is symmetric positive definite.
// The function returns the lower triangular matrix L.
void Decomposition::Cholesky(const BaseMatrix& A, BaseMatrix& L) {
    int n = A.getRows();
    assert(A.getCols() == n); // A must be square

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            // Sum L[i][k] * L[j][k]
            for (int k = 0; k < j; ++k) {
                sum += L.get(i, k) * L.get(j, k);
            }
            if (i == j) {
                // Diagonal elements: L[i][i] = sqrt(A[i][i] - sum)
                double value = A.get(i, i) - sum;
                assert(value > 0); // A must be positive definite
                L.set(i, j, std::sqrt(value));
            } else {
                // Off-diagonal elements: L[i][j] = (A[i][j] - sum) / L[j][j]
                L.set(i, j, (A.get(i, j) - sum) / L.get(j, j));
            }
        }
    }
}
