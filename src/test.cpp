#include "matrix_generator.h"
#include "base_matrix.h"
#include "EigenSolver.h"
#include <iostream>
#include <vector>

int main() {
    std::cout << "=============================" << std::endl;
    std::cout << "Sparse Matrix Library Testing" << std::endl;
    std::cout << "=============================" << std::endl;

    // ----------------------
    // Part 1: Generate a random sparse matrix and print
    // ----------------------
    MatrixGenerator mg;
    BaseMatrix* random_matrix = mg.generate_matrix("COO", 10, 10, 13);

    std::cout << "\n[Generated Random Sparse Matrix]" << std::endl;
    random_matrix->print();
    random_matrix->show_matrix();

    delete random_matrix; // Clean up random matrix

    // ----------------------
    // Part 2: Generate a symmetric positive definite (SPD) matrix
    // ----------------------
    BaseMatrix* matrix = mg.generate_spd_matrix(10);  // Generate a 10x10 SPD matrix

    std::cout << "\n[Generated SPD Matrix for EigenSolver Testing]" << std::endl;
    matrix->print();
    matrix->show_matrix();

    // ----------------------
    // Test Power Iteration
    // ----------------------
    std::cout << "\n[Testing Power Iteration]" << std::endl;
    std::pair<double, std::vector<double>> result = EigenSolver::powerIteration(*matrix);
    std::cout << "Estimated dominant eigenvalue (Power Iteration): " << result.first << std::endl;

    // ----------------------
    // Test Inverse Iteration
    // ----------------------
    std::cout << "\n[Testing Inverse Iteration]" << std::endl;
    std::pair<double, std::vector<double>> inv_result = EigenSolver::inverseIteration(*matrix);
    std::cout << "Estimated smallest eigenvalue (Inverse Iteration): " << inv_result.first << std::endl;

    // ----------------------
    // Test Lanczos Iteration
    // ----------------------
    std::vector<int> empty_row, empty_col;
    std::vector<double> empty_values;
    COO T(empty_row, empty_col, empty_values, 10, 10);
    COO Q(empty_row, empty_col, empty_values, 10, 10);

    std::cout << "\n[Testing Lanczos Iteration]" << std::endl;
    std::vector<double> lanczos_eigenvalues = EigenSolver::lanczos(*matrix, 10, T, Q);
    std::cout << "[Lanczos finished]" << std::endl;
    //T.print();

    std::cout << "[Lanczos eigenvalues]" << std::endl;
    for (double eig : lanczos_eigenvalues) {
        std::cout << eig << " ";
    }
    std::cout << std::endl;

    // ----------------------
    // Test Arnoldi Iteration
    // ----------------------
    COO H(empty_row, empty_col, empty_values, 10, 10);
    COO Q2(empty_row, empty_col, empty_values, 10, 10);

    std::cout << "\n[Testing Arnoldi Iteration]" << std::endl;
    std::vector<double> arnoldi_eigenvalues = EigenSolver::arnoldi(*matrix, 10, H, Q2);
    std::cout << "[Arnoldi finished]" << std::endl;
    //H.print();

    std::cout << "[Arnoldi eigenvalues]" << std::endl;
    for (double eig : arnoldi_eigenvalues) {
        std::cout << eig << " ";
    }
    std::cout << std::endl;


    // ----------------------
    // End
    // ----------------------
    delete matrix;

    std::cout << "\nAll tests completed." << std::endl;
    return 0;
}
