#include "matrix_generator.h"
#include "base_matrix.h"
#include "EigenSolver.h"
#include <iostream>
#include <vector>
#include <chrono> 

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
    int matrix_size = 5000; 
    BaseMatrix* matrix = mg.generate_spd_matrix(matrix_size);

    std::cout << "\n[Generated SPD Matrix for EigenSolver Testing]" << std::endl;
    //matrix->print();
    //matrix->show_matrix();

    // ----------------------
    // Test Power Iteration
    // ----------------------
    std::cout << "\n[Testing Power Iteration]" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<double, std::vector<double>> result = EigenSolver::powerIteration(*matrix);
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
    std::pair<double, std::vector<double>> inv_result = EigenSolver::inverseIteration(*matrix);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Estimated smallest eigenvalue (Inverse Iteration): " << inv_result.first << std::endl;
    std::cout << "Time taken: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;


    // ----------------------
    // Test Dense QR Iteration
    // ----------------------
    std::cout << "\n[Testing Dense QR Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> dense_eigenvalues = EigenSolver::denseQR(*matrix, matrix_size, 50);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "[Dense QR finished]" << std::endl;
    std::cout << "Time taken: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << " ms" << std::endl;

    std::cout << "[Dense QR eigenvalues (first 10)]" << std::endl;
    int print_count = std::min(10, (int)dense_eigenvalues.size());
    for (int i = 0; i < print_count; ++i) {
        std::cout << dense_eigenvalues[i] << " ";
    }
    std::cout << std::endl;



    // ----------------------
    // Test Lanczos Iteration
    // ----------------------
    std::vector<int> empty_row, empty_col;
    std::vector<double> empty_values;
    COO T(empty_row, empty_col, empty_values, matrix_size, matrix_size);
    COO Q(empty_row, empty_col, empty_values, matrix_size, matrix_size);

    std::cout << "\n[Testing Lanczos Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> lanczos_eigenvalues = EigenSolver::lanczos(*matrix, matrix_size, T, Q);
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
    COO H(empty_row, empty_col, empty_values, matrix_size, matrix_size);
    COO Q2(empty_row, empty_col, empty_values, matrix_size, matrix_size);

    std::cout << "\n[Testing Arnoldi Iteration]" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> arnoldi_eigenvalues = EigenSolver::arnoldi(*matrix, matrix_size, H, Q2);
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
    // End
    // ----------------------
    delete matrix;
    std::cout << "\nAll tests completed." << std::endl;
    return 0;
}
