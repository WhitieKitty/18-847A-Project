#include "svd.h"
#include <iostream>
#include "base_matrix.h"
#include "coo.h"
#include "matrix_generator.h"

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  extern "C" {
    #include <cblas.h>
    #include <lapacke.h>
  }
#endif

int main() {
    std::cout << "=============================" << std::endl;
    std::cout << "Sparse Matrix SVD Testing" << std::endl;
    std::cout << "=============================" << std::endl;

    MatrixGenerator mg;
    BaseMatrix* matrix = mg.generate_matrix("COO", 20, 20, 3);
    matrix->print();
    matrix->show_matrix();

    double largest_singular_value = SVD::lanczos_bidiag(*matrix, 3);
    std::cout << "Estimated Largest singular value: " << largest_singular_value << std::endl;

    std::vector<double> ori_matrx;
    for (int i = 0; i < matrix->getRows(); i++) {
        for (int j = 0; j < matrix->getCols(); j++) {
            ori_matrx.push_back(matrix->get(i, j));
        }
    }

    int m = matrix->getRows(), n = matrix->getCols(), lda = matrix->getCols();
    int ldu = m, ldvt = n;
    std::vector<double> s(std::min(m, n));
    std::vector<double> u(ldu * m);
    std::vector<double> vt(ldvt * n);

    int info;

#ifdef __APPLE__
    int lwork = -1;
    double wkopt;
    double* work;

    // Workspace query
    dgesvd_("All", "All", &m, &n,
            ori_matrx.data(), &lda,
            s.data(), u.data(), &ldu,
            vt.data(), &ldvt,
            &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    dgesvd_( "All", "All", &m, &n, ori_matrx.data(), &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info );
    std::cout << "S:";
    for (int i = 0; i < n; i++) {
        std::cout << s[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Largest singular value: " << s[0] << std::endl;
    delete matrix;
   
}
