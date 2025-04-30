#include "EigenSolver.h"
#include "coo.h"
#include "Decomposition.h"
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>




std::pair<double, std::vector<double>> EigenSolver::powerIteration(const BaseMatrix& A, int maxIters, double tol) {
    assert(A.getRows() == A.getCols());
    int n = A.getRows();
    std::vector<double> b_k(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) b_k[i] = dis(gen);
    for (int iter = 0; iter < maxIters; ++iter) {
        std::vector<double> Ab = A.multiply(b_k);
        double norm = std::sqrt(std::inner_product(Ab.begin(), Ab.end(), Ab.begin(), 0.0));

        for (int i = 0; i < n; ++i) Ab[i] /= norm;
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (Ab[i] - b_k[i]) * (Ab[i] - b_k[i]);
        diff = std::sqrt(diff);
        if (diff < tol)
        {
            std::cout << "Finish iteration in iter "<< iter<< std::endl;
            break;
        }
        b_k = Ab;
    }
    std::vector<double> Ab = A.multiply(b_k);
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

        double norm = std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));

        for (int i = 0; i < n; ++i) x[i] /= norm;
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (x[i] - b_k[i]) * (x[i] - b_k[i]);
        diff = std::sqrt(diff);
        if (diff < tol) 
        {
            std::cout << "Finish iteration in iter "<< iter<< std::endl;
            break;
        }
        b_k = x;
    }
    std::vector<double> Ab = A.multiply(b_k);
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

std::vector<double> EigenSolver::QRIteration(const BaseMatrix& A, int max_iters, double tol) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols() && "Matrix must be square.");

    COO current(A);
    
    COO Qmat({}, {}, {}, n, n);
    COO Rmat({}, {}, {}, n, n);

    for (int iter = 0; iter < max_iters; ++iter) {
        Decomposition::QR(current, Qmat, Rmat);

        // current=R@Q
        BaseMatrix* mult = Rmat.multiply(Qmat);
        COO* new_current = dynamic_cast<COO*>(mult);
        assert(new_current != nullptr);
        current = *new_current;
        delete new_current;

        double off_diag_norm = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j) off_diag_norm += current.get(i, j) * current.get(i, j);
        if (std::sqrt(off_diag_norm) < tol)
        {
            std::cout << "Finish QR iteration in iter "<< iter<< std::endl;
            break;
        }
    }

    std::vector<double> eigenvalues(n);
    for (int i = 0; i < n; ++i)
        eigenvalues[i] = current.get(i, i);
    return eigenvalues;
}

std::vector<double> EigenSolver::lanczos(const BaseMatrix& A, int k, BaseMatrix& T, BaseMatrix& Q, int max_iters) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());

    std::vector<double> q_old(n, 0.0);
    std::vector<double> q_new(n, 0.0);
    std::vector<double> w(n, 0.0);
    std::vector<double> alpha(k, 0.0);
    std::vector<double> beta(k, 0.0);

    std::mt19937 gen(777);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) q_new[i] = dis(gen);
    double norm = std::sqrt(std::inner_product(q_new.begin(), q_new.end(), q_new.begin(), 0.0));
    for (double& val : q_new) val /= norm;

    int valid_k = 0;
    for (int j = 0; j < k; ++j) {
        w = A.multiply(q_new);
        if (j > 0) {
            for (int i = 0; i < n; ++i) w[i] -= beta[j - 1] * q_old[i];
        }

        alpha[j] = std::inner_product(q_new.begin(), q_new.end(), w.begin(), 0.0);
        for (int i = 0; i < n; ++i) w[i] -= alpha[j] * q_new[i];

        beta[j] = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
        if (beta[j] < 1e-12 || std::isnan(beta[j])) {
            valid_k = j + 1;
            break;
        }

        for (int i = 0; i < n; ++i) {
            q_old[i] = q_new[i];
            q_new[i] = w[i] / beta[j];
        }
        valid_k = j + 1;
    }

    for (int i = 0; i < valid_k; ++i) {
        T.set(i, i, alpha[i]);
        if (i < valid_k - 1) {
            T.set(i, i + 1, beta[i]);
            T.set(i + 1, i, beta[i]);
        }
    }

    return QRIteration(T, max_iters, 1e-10);
}

std::vector<double> EigenSolver::arnoldi(const BaseMatrix& A, int k, BaseMatrix& H, BaseMatrix& Q, int max_iters) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());

    std::vector<std::vector<double>> q_vectors(k + 1, std::vector<double>(n, 0.0));
    std::mt19937 gen(888);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) q_vectors[0][i] = dis(gen);
    double norm = std::sqrt(std::inner_product(q_vectors[0].begin(), q_vectors[0].end(), q_vectors[0].begin(), 0.0));
    for (double& val : q_vectors[0]) val /= norm;

    int valid_k = 0;
    for (int j = 0; j < k; ++j) {
        std::vector<double> w = A.multiply(q_vectors[j]);

        for (int i = 0; i <= j; ++i) {
            double hij = std::inner_product(q_vectors[i].begin(), q_vectors[i].end(), w.begin(), 0.0);
            H.set(i, j, hij);
            for (int l = 0; l < n; ++l) w[l] -= hij * q_vectors[i][l];
        }

        double h_next = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
        
        if (h_next < 1e-12 || std::isnan(h_next)) {
            valid_k = j + 1;
            break;
        }

        if (j + 1 < k) {
            H.set(j + 1, j, h_next);
            for (int i = 0; i < n; ++i) q_vectors[j + 1][i] = w[i] / h_next;
        }
        valid_k = j + 1;
    }

    return QRIteration(H, max_iters, 1e-10);
}


