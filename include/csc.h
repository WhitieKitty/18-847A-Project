#ifndef CSC_H
#define CSC_H

#include "base_matrix.h"
#include "csc.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

/**
 * @class CSC
 * @brief Class representing a sparse matrix in Compressed Sparse Column (CSC) format.
 * 
 * Stores nonzero values column by column, with arrays for values, row indices, and column pointers.
 * It's efficient for column slicing and matrix-vector multiplication.
 */
class CSC : public BaseMatrix {
    public:
        CSC(const std::string& filename);
        CSC(const std::vector<int>& col_ptr, const std::vector<int>& row_idx, const std::vector<double>& values, int m, int n);
        ~CSC();
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
        std::vector<int> col_ptr;
        std::vector<int> row_idx;
        std::vector<double> values;
        int m; // number of rows
        int n; // number of columns
};

#endif
