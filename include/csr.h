#ifndef CSR_H
#define CSR_H
#include <vector>
#include <string> 
#include <fstream>
#include <iostream>
#include "base_matrix.h"

/**
 * @class CSR
 * @brief Class representing a sparse matrix in Compressed Sparse Row (CSR) format.
 * 
 * Stores nonzero values row by row, with arrays for values, column indices, and row start offsets. 
 * It's efficient for row slicing and matrix-vector multiplication.
 */
class CSR : public BaseMatrix {
    public:
        CSR(const std::string& filename);
        CSR(const std::vector<int>& row_ptr, const std::vector<int>& col_idx, const std::vector<double>& values, int m, int n);
        ~CSR();
        void print() const;
        void show_matrix() const;
        BaseMatrix* multiply(const BaseMatrix& other) const;
        std::vector<double> multiply(const std::vector<double>& other) const;
        BaseMatrix* add(const BaseMatrix& other) const;
        BaseMatrix* subtract(const BaseMatrix& other) const;
        BaseMatrix* transpose() const;
        int getRows() const;
        int getCols() const;
        void set(int i, int j, double value);
        double get(int i, int j) const;
    private:
        std::vector<int> row_ptr;    // Array of row pointers (start index of each row)
        std::vector<int> col_idx;    // Array of column indices for non-zero elements
        std::vector<double> values;  // Array of non-zero values
        int m; // number of rows
        int n; // number of columns
};

#endif
