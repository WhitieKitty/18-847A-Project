#include "matrix_generator.h"
#include "base_matrix.h"
#include "EigenSolver.h"
#include "Decomposition.h"
#include <iostream>
#include <vector>
#include <chrono> 

double frobenius_norm(BaseMatrix* mat) {
    double frobenius_norm = 0.0;
    for (int i = 0; i < mat->getRows(); ++i) {
        for (int j = 0; j < mat->getCols(); ++j) {
            double value = mat->get(i, j);
            frobenius_norm += value * value;
        }
    }
    return std::sqrt(frobenius_norm);
}

int main() {
    std::cout << "=============================" << std::endl;
    std::cout << "Sparse Matrix Library Testing" << std::endl;
    std::cout << "=============================" << std::endl;

    

    // ----------------------
    // Part 1: Generate a random sparse matrix and print
    // ----------------------
    MatrixGenerator mg;
    BaseMatrix* random_matrix = mg.generate_matrix("COO", 10, 10, 3);

    std::cout << "\n[Generated Random Sparse Matrix]" << std::endl;
    random_matrix->print();
    random_matrix->show_matrix();

    delete random_matrix; 

    // ----------------------
    // Part 2: Generate a symmetric positive definite (SPD) matrix
    // ----------------------
    int matrix_size = 50; // the size of spd matrix
    int max_iter = 50; // the maximum iterations for solving eigen values
    int num_eigenvalues=10; // the number of eigenvalues to be calculated for Lanczos and Arnoldi 
    
    BaseMatrix* matrix = mg.generate_spd_matrix("COO", matrix_size);

    std::cout << "\n[Generated SPD Matrix for EigenSolver Testing]" << std::endl;
    //matrix->print();
    //matrix->show_matrix();

    // ----------------------
    // Test Power Iteration
    // ----------------------
    std::cout << "\n[Testing Power Iteration]" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<double, std::vector<double>> result = EigenSolver::powerIteration(*matrix,max_iter,1e-10);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Estimated dominant eigenvalue (Power Iteration): " << result.first << std::endl;
    std::cout << "Time taken: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;

    // ----------------------
    // Test Inverse Iteration
    // ----------------------
    std::cout << "\n[Testing Inverse Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::pair<double, std::vector<double>> inv_result = EigenSolver::inverseIteration(*matrix,max_iter,1e-10);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Estimated smallest eigenvalue (Inverse Iteration): " << inv_result.first << std::endl;
    std::cout << "Time taken: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;


    // ----------------------
    // Test QR Iteration
    // ----------------------
    std::cout << "\n[Testing QR Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> qr_eigenvalues = EigenSolver::QRIteration(*matrix, max_iter, 1e-10);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[QR finished]" << std::endl;
    std::cout << "Time taken: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << " ms" << std::endl;

    std::cout << "[QR eigenvalues (first 10)]" << std::endl;
    int print_count = std::min(10, (int)qr_eigenvalues.size());
    for (int i = 0; i < print_count; ++i) {
        std::cout << qr_eigenvalues[i] << " ";
    }
    std::cout << std::endl;



    // ----------------------
    // Test Lanczos Iteration
    // ----------------------
    std::vector<int> empty_row, empty_col;
    std::vector<double> empty_values;
    COO T(empty_row, empty_col, empty_values, num_eigenvalues, num_eigenvalues);
    COO Q(empty_row, empty_col, empty_values, num_eigenvalues, num_eigenvalues);

    std::cout << "\n[Testing Lanczos Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> lanczos_eigenvalues = EigenSolver::lanczos(*matrix, num_eigenvalues, T, Q,max_iter);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Lanczos finished]" << std::endl;
    std::cout << "Time taken: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;


    std::cout << "[Lanczos eigenvalues (first 10)]" << std::endl;
    print_count = std::min(10, (int)lanczos_eigenvalues.size());
    for (int i = 0; i < print_count; ++i) {
        std::cout << lanczos_eigenvalues[i] << " ";
    }
    std::cout << std::endl;
              

    // ----------------------
    // Test Arnoldi Iteration
    // ----------------------
    COO H(empty_row, empty_col, empty_values, num_eigenvalues, num_eigenvalues);
    COO Q2(empty_row, empty_col, empty_values, num_eigenvalues, num_eigenvalues);

    std::cout << "\n[Testing Arnoldi Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> arnoldi_eigenvalues = EigenSolver::arnoldi(*matrix, num_eigenvalues, H, Q2,max_iter);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Arnoldi finished]" << std::endl;
    std::cout << "Time taken: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;


    std::cout << "[Arnoldi eigenvalues (first 10)]" << std::endl;
    print_count = std::min(10, (int)arnoldi_eigenvalues.size());
    for (int i = 0; i < print_count; ++i) {
        std::cout << arnoldi_eigenvalues[i] << " ";
    }
    std::cout << std::endl;

    // ----------------------
    // Part 3: Test Decomposition
    // ----------------------
    std::cout << "\n=========================================" << std::endl;
    std::cout << "        Decomposition Testing" << std::endl;
    std::cout << "=========================================" << std::endl;

    // ----------------------
    // Test LU Decomposition
    // ----------------------
    std::cout << "\n[Generated Random Sparse Matrix in COO Format for LU Decomposition]" << std::endl;
    matrix_size = 10;
    BaseMatrix* lu_matrix = mg.generate_spd_matrix("COO", matrix_size);
    lu_matrix->print();

    BaseMatrix* L = mg.generate_matrix("COO", matrix_size, matrix_size, 0);  // empty matrix for L
    BaseMatrix* U = mg.generate_matrix("COO", matrix_size, matrix_size, 0);  // empty matrix for U

    start = std::chrono::high_resolution_clock::now();
    Decomposition::LU(*lu_matrix, *L, *U);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "[L Matrix]" << std::endl;
    L->show_matrix();
    std::cout << "[U Matrix]" << std::endl;
    U->show_matrix();

    auto lu_product = L->multiply(*U);

    // std::cout << "[L * U Product]" << std::endl;
    // lu_product->print();

    auto difference = lu_product->subtract(*lu_matrix);
    double fb_norm = frobenius_norm(difference);
    std::cout << "[Frobenius Norm of Difference (L * U - Original Matrix)]" << std::endl;
    std::cout << "Frobenius Norm: " << fb_norm << std::endl;

    std::cout << "LU Decomposition Time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
              << " ns" << std::endl;

    delete lu_matrix;
    delete L;
    delete U;

    // ----------------------
    // Test QR Decomposition
    // ----------------------
    std::cout << "\n[Generated Random Sparse Matrix in CSR Format for QR Decomposition]" << std::endl;
    BaseMatrix* qr_matrix = mg.generate_spd_matrix("CSR", matrix_size);
    qr_matrix->print();

    BaseMatrix* Q_mat = mg.generate_matrix("CSR", matrix_size, matrix_size, 0);  // empty matrix for Q
    BaseMatrix* R = mg.generate_matrix("CSR", matrix_size, matrix_size, 0);  // empty matrix for R

    start = std::chrono::high_resolution_clock::now();
    Decomposition::QR(*qr_matrix, *Q_mat, *R);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "[Q Matrix]" << std::endl;
    Q_mat->show_matrix();
    std::cout << "[R Matrix]" << std::endl;
    R->show_matrix();
    auto qr_product = Q_mat->multiply(*R);
    auto difference_qr = qr_product->subtract(*qr_matrix);
    double frobenius_norm_qr = frobenius_norm(difference_qr);
    std::cout << "Frobenius Norm of Difference (Q * R - Original Matrix)" << std::endl;
    std::cout << "Frobenius Norm: " << frobenius_norm_qr << std::endl;


    std::cout << "QR Decomposition Time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
              << " ns" << std::endl;

    delete qr_matrix;
    delete Q_mat;
    delete R;

    // ----------------------
    // Test Cholesky Decomposition
    // ----------------------
    std::cout << "\n[Generated Random Sparse Matrix in CSC Format for Cholesky Decomposition]" << std::endl;
    BaseMatrix* cholesky_matrix = mg.generate_spd_matrix("CSC", matrix_size);
    cholesky_matrix->print();

    BaseMatrix* chol_L = mg.generate_matrix("CSC", matrix_size, matrix_size, 0);  // empty matrix for L

    start = std::chrono::high_resolution_clock::now();
    Decomposition::Cholesky(*cholesky_matrix, *chol_L);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "[Cholesky L Matrix]" << std::endl;
    chol_L->show_matrix();
    auto chol_L_trans = chol_L->transpose();
    auto chol_product = chol_L->multiply(*chol_L_trans);
    auto diff_chol = chol_product->subtract(*cholesky_matrix);
    double frobenius_norm_chol = frobenius_norm(diff_chol);
    std::cout << "[Frobenius Norm of Difference (L * L^T - Original Matrix)]" << std::endl;
    std::cout << "Frobenius Norm: " << frobenius_norm_chol << std::endl;

    std::cout << "Cholesky Decomposition Time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
              << " ns" << std::endl;

    delete cholesky_matrix;
    delete chol_L;

    std::cout << "\nAll decomposition tests completed." << std::endl;


    // ----------------------
    // End
    // ----------------------
    delete matrix;
    std::cout << "\nAll tests completed." << std::endl;
    return 0;
}
