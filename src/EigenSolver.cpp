#include "EigenSolver.h"
#include "coo.h"
#include "Decomposition.h"
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>



// Power Iteration Method
// This method estimates the largest (dominant) eigenvalue and its corresponding eigenvector of a square matrix A.
// Steps:
// (1) Initialize a random starting vector b_0
// (2) Iteratively update b_{k+1} = A * b_k and normalize it
//     - Stop the iteration when the update converges (i.e., ||b_{k+1} - b_k|| < tol)
// (3) Estimate the dominant eigenvalue using the Rayleigh quotient: λ ≈ (b_k^T A b_k) / (b_k^T b_k)
std::pair<double, std::vector<double>> EigenSolver::powerIteration(const BaseMatrix& A, int maxIters, double tol) {
    // Assumes A is square: n x n
    assert(A.getRows() == A.getCols());
    int n = A.getRows();

    // Step 1: initialize a random vector b_0
    std::vector<double> b_k(n);  // b_k is the vector for current iteration
    std::mt19937 gen(42); // random number generator with fixed seed for reproducibility
    std::uniform_real_distribution<> dis(0.0, 1.0); // uniform distribution in [0, 1]
    for (int i = 0; i < n; ++i) b_k[i] = dis(gen);  // fill b_k with random values

    // Step 2: iterate b_{k+1} = A b_k, then normalize and check convergence
    for (int iter = 0; iter < maxIters; ++iter) 
    {
        std::vector<double> Ab = A.multiply(b_k); // compute A * b_k
        double norm = std::sqrt(std::inner_product(Ab.begin(), Ab.end(), Ab.begin(), 0.0));  // ||A b_k||

        for (int i = 0; i < n; ++i) Ab[i] /= norm; // normalize Ab to unit length (b_{k+1})
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (Ab[i] - b_k[i]) * (Ab[i] - b_k[i]);
        diff = std::sqrt(diff); // compute ||b_{k+1} - b_k||
        if (diff < tol)
        {
            std::cout << "Finish iteration in iter "<< iter<< std::endl;
            break;  // convergence criterion met
        }
        b_k = Ab;  // update vector for next iteration
    }

    // Step 3: Rayleigh quotient: lambda ≈ (b_k^T A b_k) / (b_k^T b_k)
    std::vector<double> Ab = A.multiply(b_k);   // Ab = A * b_k for final Rayleigh quotient
    double numerator = 0.0, denominator = 0.0;
    for (int i = 0; i < n; ++i) 
    {
        numerator += b_k[i] * Ab[i];    // b_k^T A b_k
        denominator += b_k[i] * b_k[i];   // b_k^T b_k
    }
    double eigenvalue = numerator / denominator;   // estimated dominant eigenvalue
    return {eigenvalue, b_k};   // return dominant eigenvalue and eigenvector
}


// Inverse Iteration Method
// This method estimates the eigenvalue closest to zero (or to a shift u if modified).
// Steps:
// (1) Initialize a random starting vector b_0
// (2) Perform LU decomposition: A = LU
// (3) Iterate: Solve LU x = b_k (i.e., x ≈ A⁻¹ b_k), normalize x, and check convergence
// (4) Estimate eigenvalue using Rayleigh quotient: λ ≈ (b_k^T A b_k) / (b_k^T b_k)

std::pair<double, std::vector<double>> EigenSolver::inverseIteration(const BaseMatrix& A, int maxIters, double tol) {
    assert(A.getRows() == A.getCols());
    int n = A.getRows();

    // Step 1: initialize a random vector b_0
    std::vector<double> b_k(n); // b_k is the vector for current iteration
    std::mt19937 gen(123); // random number generator with fixed seed for reproducibility
    std::uniform_real_distribution<> dis(0.0, 1.0); // uniform distribution in [0, 1]
    for (int i = 0; i < n; ++i) b_k[i] = dis(gen); // fill b_k with random values

    // Step 2: LU decomposition A = LU
    BaseMatrix* L = new COO({}, {}, {}, n, n);
    BaseMatrix* U = new COO({}, {}, {}, n, n);
    Decomposition::LU(A, *L, *U);

    // Step 3: inverse iteration,iterate b_{k+1} = A^{-1} b_k, then normalize and check convergence
    for (int iter = 0; iter < maxIters; ++iter) 
    {
        // Solve the linear system LU x = b_k, which approximates x ≈ A^{-1} b_k
        std::vector<double> x = Decomposition::solveLU(*L, *U, b_k);

        // Compute the Euclidean norm ||x|| to normalize the solution vector
        double norm = std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
        for (int i = 0; i < n; ++i) x[i] /= norm; // Normalize the vector

        // Compute difference ||x - b_k|| to check for convergence
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += (x[i] - b_k[i]) * (x[i] - b_k[i]);
        diff = std::sqrt(diff);

        // If the vector is sufficiently close to the previous iteration, stop
        if (diff < tol) {
            std::cout << "Finish iteration in iter " << iter << std::endl;
            break;
        }

        // Update b_k to the new normalized vector x
        b_k = x;
    }

    // Step 4: estimate eigenvalue using Rayleigh quotient λ ≈ (b_k^T A b_k) / (b_k^T b_k)
    std::vector<double> Ab = A.multiply(b_k); // Compute A * b_k for numerator
    double numerator = 0.0, denominator = 0.0;
    for (int i = 0; i < n; ++i) 
    {
        numerator += b_k[i] * Ab[i];      // Dot product: b_k^T A b_k
        denominator += b_k[i] * b_k[i];   // Dot product: b_k^T b_k (squared norm)
    }

    // Final estimated eigenvalue closest to 0 (or shift)
    double eigenvalue = numerator / denominator;
    delete L;
    delete U;
    return {eigenvalue, b_k};
}


// QR Iteration Method
// This method estimates all eigenvalues of a square matrix A using repeated QR decomposition.
//
// Steps:
// (1) Let current matrix A_0 = A
// (2) Repeat the QR iteration:
//     - Compute QR decomposition of A_k → A_k = Q_k R_k
//     - Update A_{k+1} = R_k Q_k
//     - This preserves eigenvalues and shifts matrix toward upper triangular form
//     - Stop when off-diagonal entries are sufficiently small
// (3) After convergence, eigenvalues are approximated by diagonal elements of A_k
//
// Why it works:
// - A_k+1 = Q_kT A_k Q_k is a similarity transform → preserves eigenvalues
// - Converges to diagonal matrix for symmetric matrices


std::vector<double> EigenSolver::QRIteration(const BaseMatrix& A, int max_iters, double tol) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols() && "Matrix must be square.");

    COO current(A); // A_0 = A
    COO Qmat({}, {}, {}, n, n); // Holds Q
    COO Rmat({}, {}, {}, n, n); // Holds R
    //Repeat the QR iteration:
    for (int iter = 0; iter < max_iters; ++iter) 
    {
        Decomposition::QR(current, Qmat, Rmat); // Step 1: QR decomposition of A_k
        BaseMatrix* mult = Rmat.multiply(Qmat); // Step 2: form A_{k+1} = R_k Q_k
        COO* new_current = dynamic_cast<COO*>(mult);
        assert(new_current != nullptr);
        current = *new_current; // Replace current with A_{k+1}
        delete new_current;

        // Step 3: check convergence by computing norm of off-diagonal entries
        double off_diag_norm = 0.0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (i != j) off_diag_norm += current.get(i, j) * current.get(i, j);
            }    
        }    
        if (std::sqrt(off_diag_norm) < tol) 
        {
            std::cout << "Finish QR iteration in iter " << iter << std::endl;
            break;
        }
    }

    // Step 4: return the diagonal entries as eigenvalues
    std::vector<double> eigenvalues(n);
    for (int i = 0; i < n; ++i)
        eigenvalues[i] = current.get(i, i);
    return eigenvalues;
}



// Lanczos Iteration Method
// This method approximates the k largest (or smallest) eigenvalues of a symmetric real matrix A ∈ ℝⁿˣⁿ.
// It is particularly efficient for large sparse matrices.
// 
// Core idea:
// - Construct an orthonormal basis Qₖ = [q₀, ..., qₖ₋₁] for the Krylov subspace:
//     Kₖ(A, q₀) = span{q₀, Aq₀, A²q₀, ..., A^{k-1}q₀}
// - Project A onto this subspace to form a small tridiagonal matrix Tₖ = Qₖᵀ A Qₖ
// - A ≈ Qₖ Tₖ Qₖᵀ ⇒ eigenvalues of Tₖ approximate those of A
// - Use QR iteration on Tₖ to estimate eigenvalues
// 
// The recurrence used to maintain orthogonality and tridiagonal structure is:
//     A qⱼ = βⱼ₋₁ qⱼ₋₁ + αⱼ qⱼ + βⱼ qⱼ₊₁
//     where:
//     - αⱼ = qⱼᵀ A qⱼ (diagonal of Tₖ)
//     - βⱼ = ||w|| (off-diagonal of Tₖ)
//     - w is the orthogonalized residual

std::vector<double> EigenSolver::lanczos(const BaseMatrix& A, int k, BaseMatrix& T, BaseMatrix& Q, int max_iters) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());

    std::vector<double> q_old(n, 0.0); // q_{j-1}, initialized to 0
    std::vector<double> q_new(n, 0.0); // q_j, the current basis vector
    std::vector<double> w(n, 0.0);     // residual vector
    std::vector<double> alpha(k, 0.0); // diagonal entries α_j
    std::vector<double> beta(k, 0.0);  // off-diagonal entries β_j

    // Step 1: Initialize q₀ to a random normalized vector
    std::mt19937 gen(777); // random number generator with fixed seed for reproducibility
    std::uniform_real_distribution<> dis(0.0, 1.0); // Uniform distribution in range [0, 1] to fill vector with random values

    // Generate initial random vector q_new with values in [0,1]
    for (int i = 0; i < n; ++i) q_new[i] = dis(gen);

    // Normalize the initial vector: q_new ← q_new / ||q_new||
    // Ensures first Lanczos vector q₀ has unit norm, forming orthonormal basis Qₖ
    double norm = std::sqrt(std::inner_product(q_new.begin(), q_new.end(), q_new.begin(), 0.0));
    for (double& val : q_new) val /= norm;

    int valid_k = 0; // Number of successfully constructed Lanczos vectors
    for (int j = 0; j < k; ++j) 
    {
        // Step 2: Compute w = A * q_j
        w = A.multiply(q_new);

        // Step 3: Subtract previous β term * q_{j-1}
        if (j > 0) 
        {
            for (int i = 0; i < n; ++i) w[i] -= beta[j - 1] * q_old[i];
        }

        // Step 4: Compute α_j = q_jᵀ w
        alpha[j] = std::inner_product(q_new.begin(), q_new.end(), w.begin(), 0.0);

        // Step 5: Orthogonalize: w ← w - α_j * q_j
        for (int i = 0; i < n; ++i) w[i] -= alpha[j] * q_new[i];

        // Step 6: Compute β_j = ||w||
        beta[j] = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));

        // Step 7: Stop if β is too small or NaN (loss of orthogonality)
        if (beta[j] < 1e-12 || std::isnan(beta[j])) 
        {
            valid_k = j + 1;
            break;
        }

        // Step 8: Update q_{j-1} and normalize new q_{j+1}
        for (int i = 0; i < n; ++i) 
        {
            q_old[i] = q_new[i];
            q_new[i] = w[i] / beta[j];
        }
        valid_k = j + 1;
    }

    // Step 9: Fill the symmetric tridiagonal matrix Tₖ
    for (int i = 0; i < valid_k; ++i) 
    {
        T.set(i, i, alpha[i]);
        if (i < valid_k - 1) 
        {
            T.set(i, i + 1, beta[i]);
            T.set(i + 1, i, beta[i]);
        }
    }

    // Step 10: Use QR Iteration on T to approximate eigenvalues of A
    return QRIteration(T, max_iters, 1e-10);
}


// Arnoldi Iteration Method
// This method approximates a few eigenvalues of a general (non-symmetric) square matrix A ∈ ℝⁿˣⁿ.
// It builds an orthonormal basis for the Krylov subspace and projects A into a smaller upper Hessenberg matrix Hₖ.
//
// Core idea:
// - Construct orthonormal basis Qₖ = [q₀, ..., qₖ₋₁] for Krylov subspace:
//     Kₖ(A, q₀) = span{q₀, Aq₀, A²q₀, ..., A^{k-1}q₀}
// - Compute Hₖ = Qₖᵀ A Qₖ → an upper Hessenberg matrix
// - A ≈ Qₖ Hₖ Qₖᵀ ⇒ eig(Hₖ) ≈ partial eig(A)
//
// Main loop implements the recurrence:
//     A qⱼ = ∑₀ⱼ hᵢⱼ qᵢ + h_{j+1,j} q_{j+1}
//     - hᵢⱼ = qᵢᵀ A qⱼ → projection coefficient
//     - Residual is orthogonalized against all previous qᵢ

std::vector<double> EigenSolver::arnoldi(const BaseMatrix& A, int k, BaseMatrix& H, BaseMatrix& Q, int max_iters) {
    int n = A.getRows();
    assert(A.getRows() == A.getCols());

    // Step 1: Prepare orthonormal basis storage: q₀,...,qₖ (columns of Q)
    std::vector<std::vector<double>> q_vectors(k + 1, std::vector<double>(n, 0.0));

    // Step 2: Initialize q₀ with random unit norm vector
    std::mt19937 gen(888); // random number generator with fixed seed for reproducibility
    std::uniform_real_distribution<> dis(0.0, 1.0); // uniform distribution [0,1]
    for (int i = 0; i < n; ++i) q_vectors[0][i] = dis(gen);
    double norm = std::sqrt(std::inner_product(q_vectors[0].begin(), q_vectors[0].end(), q_vectors[0].begin(), 0.0));
    for (double& val : q_vectors[0]) val /= norm; // normalize q₀

    int valid_k = 0; // number of Arnoldi steps performed
    for (int j = 0; j < k; ++j) 
    {
        std::vector<double> w = A.multiply(q_vectors[j]); // w = A * qⱼ

        // Orthogonalize w against all previous qᵢ
        for (int i = 0; i <= j; ++i) 
        {
            double hij = std::inner_product(q_vectors[i].begin(), q_vectors[i].end(), w.begin(), 0.0); // hᵢⱼ = qᵢᵀ w
            H.set(i, j, hij); // store in Hessenberg matrix
            for (int l = 0; l < n; ++l) w[l] -= hij * q_vectors[i][l]; // w -= hᵢⱼ * qᵢ
        }

        // Compute h_{j+1, j} = ||w||
        double h_next = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));

        // Stop if breakdown or converged (null residual)
        if (h_next < 1e-12 || std::isnan(h_next)) 
        {
            valid_k = j + 1;
            break;
        }

        // Normalize q_{j+1} = w / h_{j+1, j}
        if (j + 1 < k) 
        {
            H.set(j + 1, j, h_next);
            for (int i = 0; i < n; ++i) q_vectors[j + 1][i] = w[i] / h_next;
        }
        valid_k = j + 1;
    }

    // Step 3: Compute eigenvalues of Hₖ using QR iteration
    return QRIteration(H, max_iters, 1e-10);
}
