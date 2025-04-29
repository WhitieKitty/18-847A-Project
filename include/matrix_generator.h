#ifndef MATRIX_GENERATOR_H
#define MATRIX_GENERATOR_H
#include <vector>
#include <string>
#include "base_matrix.h"
#include "csr.h"
#include "csc.h"
#include "coo.h"
#include <random>
#include <algorithm>

class MatrixGenerator {
    public:
        MatrixGenerator();
        ~MatrixGenerator();
        BaseMatrix* generate_matrix(const std::string& format, int m, int n, int density);
        BaseMatrix* generate_spd_matrix(int n);

        private:
};

#endif
