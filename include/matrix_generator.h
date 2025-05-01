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

/**
 * @class MatrixGenerator
 * @brief A utility class for generating random sparse matrices in various formats.
 * 
 */
class MatrixGenerator {
    public:
        MatrixGenerator();
        ~MatrixGenerator();

        /**
         * @brief Generates a random sparse matrix in the specified format.
         * 
         * @param format The format of the matrix (e.g., "CSR", "CSC", "COO").
         * @param m Number of rows.
         * @param n Number of columns.
         * @param density Percentage of non-zero elements (0-100).
         * 
         * @return A pointer to the generated sparse matrix in the specified format.
         */
        BaseMatrix* generate_matrix(const std::string& format, int m, int n, int density);

        /**
         * @brief Generates a symmetric positive definite (SPD) matrix in the specified format.
         * 
         * @param format The format of the matrix (e.g., "CSR", "CSC", "COO").
         * @param n Size of the matrix (n x n).
         * 
         * @return A pointer to the generated SPD matrix in the specified format.
         */
        BaseMatrix* generate_spd_matrix(const std::string& format, int n);

        private:
};

#endif
