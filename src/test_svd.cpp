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
    work = (double*)malloc(lwork * sizeof(double));

    // Actual computation
    dgesvd_("All", "All", &m, &n,
            ori_matrx.data(), &lda,
            s.data(), u.data(), &ldu,
            vt.data(), &ldvt,
            work, &lwork,
            &info);

    free(work);
#else
    std::vector<double> superb(std::min(m, n) - 1);
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A',
                          m, n, ori_matrx.data(), n,
                          s.data(), u.data(), ldu,
                          vt.data(), ldvt, superb.data());
#endif

    if (info != 0) {
        std::cerr << "SVD failed, info = " << info << std::endl;
    } else {
        std::cout << "S: ";
        for (auto val : s) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "Largest singular value: " << s[0] << std::endl;
    }

    delete matrix;
}
