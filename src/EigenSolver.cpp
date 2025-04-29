#include "EigenSolver.h"
#include "coo.h"
#include "Decomposition.h"
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>

std::vector<double> sparse_matvec_product(const BaseMatrix& A, const std::vector<double>& x) {
    int rows = A.getRows();
    int cols = A.getCols();
    assert((int)x.size() == cols && "Dimension mismatch in sparse_matvec_product.");
    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double aij = A.get(i, j);
            if (std::abs(aij) > 1e-12) {
                result[i] += aij * x[j];
            }
        }
    }
    return result;
}

std::pair<double, std::vector<double>> EigenSolver::powerIteration(const BaseMatrix& A, int maxIters, double tol) {
    assert(A.getRows() == A.getCols());
    int n = A.getRows();
    std::vector<double> b_k(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) b_k[i] = dis(gen);
    for (int iter = 0; iter < maxIters; ++iter) {
        std::vector<double> Ab = sparse_matvec_product(A, b_k);
        double norm = 0.0;
        for (double val : Ab) norm += val * val;
        norm = std::sqrt(norm);
        for (int i = 0; i < n; ++i) Ab[i] /= norm;
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (Ab[i] - b_k[i]) * (Ab[i] - b_k[i]);
        diff = std::sqrt(diff);
        if (diff < tol) break;
        b_k = Ab;
    }
    std::vector<double> Ab = sparse_matvec_product(A, b_k);
    double numerator = 0.0, denominator = 0.0;
    for (int i = 0; i < n; ++i) {
        numerator += b_k[i] * Ab[i];
        denominator += b_k[i] * b_k[i];
    }
    double eigenvalue = numerator / denominator;
    return {eigenvalue, b_k};
}

std::pair<double, std::vector<double>> EigenSolver::inverseIteration(const BaseMatrix& A, int maxIters, double tol) {
    assert(A.getRows() == A.getCols());
    int n = A.getRows();
    std::vector<double> b_k(n);
    std::mt19937 gen(123);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) b_k[i] = dis(gen);
    BaseMatrix* L = new COO({}, {}, {}, n, n);
    BaseMatrix* U = new COO({}, {}, {}, n, n);
    Decomposition::LU(A, *L, *U);
    for (int iter = 0; iter < maxIters; ++iter) {
        std::vector<double> x = Decomposition::solveLU(*L, *U, b_k);
        double norm = 0.0;
        for (double val : x) norm += val * val;
        norm = std::sqrt(norm);
        for (int i = 0; i < n; ++i) x[i] /= norm;
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (x[i] - b_k[i]) * (x[i] - b_k[i]);
        diff = std::sqrt(diff);
        if (diff < tol) break;
        b_k = x;
    }
    std::vector<double> Ab = sparse_matvec_product(A, b_k);
    double numerator = 0.0, denominator = 0.0;
    for (int i = 0; i < n; ++i) {
        numerator += b_k[i] * Ab[i];
        denominator += b_k[i] * b_k[i];
    }
    double eigenvalue = numerator / denominator;
    delete L;
    delete U;
    return {eigenvalue, b_k};
}

std::vector<double> EigenSolver::lanczos(const BaseMatrix& A, int k, BaseMatrix& T, BaseMatrix& Q) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());
    std::vector<std::vector<double>> q_vectors(k + 1, std::vector<double>(n, 0.0));
    std::vector<double> alpha(k, 0.0);
    std::vector<double> beta(k, 0.0);
    std::mt19937 gen(777);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) q_vectors[0][i] = dis(gen);
    double norm = 0.0;
    for (double val : q_vectors[0]) norm += val * val;
    norm = std::sqrt(norm);
    for (double& val : q_vectors[0]) val /= norm;
    for (int j = 0; j < k; ++j) {
        std::vector<double> w = sparse_matvec_product(A, q_vectors[j]);
        if (j > 0) for (int i = 0; i < n; ++i) w[i] -= beta[j-1] * q_vectors[j-1][i];
        alpha[j] = 0.0;
        for (int i = 0; i < n; ++i) alpha[j] += w[i] * q_vectors[j][i];
        for (int i = 0; i < n; ++i) w[i] -= alpha[j] * q_vectors[j][i];
        beta[j] = 0.0;
        for (int i = 0; i < n; ++i) beta[j] += w[i] * w[i];
        beta[j] = std::sqrt(beta[j]);
        if (j + 1 < k) {
            if (beta[j] > 1e-12) {
                for (int i = 0; i < n; ++i) q_vectors[j+1][i] = w[i] / beta[j];
            } else {
                for (int i = 0; i < n; ++i) q_vectors[j+1][i] = dis(gen);
                norm = 0.0;
                for (double val : q_vectors[j+1]) norm += val * val;
                norm = std::sqrt(norm);
                for (double& val : q_vectors[j+1]) val /= norm;
            }
        }
    }
    for (int i = 0; i < k; ++i) {
        T.set(i, i, alpha[i]);
        if (i < k-1) {
            T.set(i, i+1, beta[i]);
            T.set(i+1, i, beta[i]);
        }
    }
    for (int j = 0; j < k; ++j) for (int i = 0; i < n; ++i) Q.set(i, j, q_vectors[j][i]);
    std::vector<double> eigenvalues(k);
    COO current_T = static_cast<const COO&>(T);
    COO Qmat({}, {}, {}, k, k);
    COO Rmat({}, {}, {}, k, k);
    for (int iter = 0; iter < 50; ++iter) {
        Decomposition::QR(current_T, Qmat, Rmat);
        COO new_T({}, {}, {}, k, k);
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                for (int l = 0; l < k; ++l)
                    new_T.set(i, j, new_T.get(i, j) + Rmat.get(i, l) * Qmat.get(l, j));
        current_T = new_T;
    }
    for (int i = 0; i < k; ++i) eigenvalues[i] = current_T.get(i, i);
    return eigenvalues;
}

std::vector<double> EigenSolver::arnoldi(const BaseMatrix& A, int k, BaseMatrix& H, BaseMatrix& Q) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());
    std::vector<std::vector<double>> q_vectors(k + 1, std::vector<double>(n, 0.0));
    std::mt19937 gen(888);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) q_vectors[0][i] = dis(gen);
    double norm = 0.0;
    for (double val : q_vectors[0]) norm += val * val;
    norm = std::sqrt(norm);
    for (double& val : q_vectors[0]) val /= norm;
    for (int j = 0; j < k; ++j) {
        std::vector<double> w = sparse_matvec_product(A, q_vectors[j]);
        for (int i = 0; i <= j; ++i) {
            double hij = 0.0;
            for (int l = 0; l < n; ++l) hij += q_vectors[i][l] * w[l];
            H.set(i, j, hij);
            for (int l = 0; l < n; ++l) w[l] -= hij * q_vectors[i][l];
        }
        double h_next = 0.0;
        for (double val : w) h_next += val * val;
        h_next = std::sqrt(h_next);
        if (h_next < 1e-12) break;
        if (j + 1 < k) {
            H.set(j+1, j, h_next);
            for (int i = 0; i < n; ++i) q_vectors[j+1][i] = w[i] / h_next;
        }
    }
    for (int j = 0; j < k; ++j) for (int i = 0; i < n; ++i) Q.set(i, j, q_vectors[j][i]);
    std::vector<double> eigenvalues(k);
    COO current_H = static_cast<const COO&>(H);
    COO Qmat({}, {}, {}, k, k);
    COO Rmat({}, {}, {}, k, k);
    for (int iter = 0; iter < 50; ++iter) {
        Decomposition::QR(current_H, Qmat, Rmat);
        COO new_H({}, {}, {}, k, k);
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                for (int l = 0; l < k; ++l)
                    new_H.set(i, j, new_H.get(i, j) + Rmat.get(i, l) * Qmat.get(l, j));
        current_H = new_H;
    }
    for (int i = 0; i < k; ++i) eigenvalues[i] = current_H.get(i, i);
    return eigenvalues;
}