#include "svd.h"
#include <iostream>
#include "base_matrix.h"
#include "coo.h"
#include "matrix_generator.h"
#include <Accelerate/Accelerate.h>

int main(int argc, char* argv[]) {
    std::cout << "=============================" << std::endl;
    std::cout << "Sparse Matrix SVD Testing" << std::endl;
    std::cout << "=============================" << std::endl;

    MatrixGenerator mg;
    int m_size, n_size, density;
    // get argv[1], argv[2], argv[3]
    if (argc != 4) {
        m_size = 20;
        n_size = 20;
        density = 3;
    } else {
        m_size = atoi(argv[1]);
        n_size = atoi(argv[2]);
        density = atoi(argv[3]);
    }
    BaseMatrix* matrix = mg.generate_matrix("COO", m_size, n_size, density);
    // matrix->print();
    // matrix->show_matrix();
    std::cout << "=============================" << std::endl;
    std::cout << "Lanczos Bidiagonalization" << std::endl;
    std::cout << "=============================" << std::endl;
    // time the lanczos_bidiag
    auto start = std::chrono::high_resolution_clock::now();
    double largest_singular_value = SVD::lanczos_bidiag(*matrix, 3);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Estimated Largest singular value: " << largest_singular_value << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0 << " ms" << std::endl;

    std::vector<double> ori_matrx;
    for (int i = 0; i < matrix->getRows(); i++) {
        for (int j = 0; j < matrix->getCols(); j++) {
            ori_matrx.push_back(matrix->get(i, j));
        }
    }
    int info;
    int m = matrix->getRows(), n = matrix->getCols(), lda = matrix->getRows();
    int ldu = matrix->getRows(), ldvt = matrix->getCols();
    double s[n], u[ldu*m], vt[ldvt*n];
    int lwork = -1;
    double wkopt;
    double* work;
    start = std::chrono::high_resolution_clock::now();
    dgesvd_( "All", "All", &m, &n, ori_matrx.data(), &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );  
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    dgesvd_( "All", "All", &m, &n, ori_matrx.data(), &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info );
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0 << " ms" << std::endl;
    // std::cout << "S:";
    // for (int i = 0; i < n; i++) {
    //     std::cout << s[i] << " ";
    // }
    // std::cout << std::endl;
    std::cout << "Largest singular value: " << s[0] << std::endl;
    delete matrix;
   
}