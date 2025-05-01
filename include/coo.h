#ifndef COO_H
#define COO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "base_matrix.h"

/**
 * @class COO
 * @brief Class representing a sparse matrix in Coordinate List (COO) format.
 * 
 * COO format stores a sparse matrix as three arrays: row indices, column indices, and values.
 */
class COO : public BaseMatrix {
    public:
        COO(const std::string& filename);
        COO(const std::vector<int>& row_idx, const std::vector<int>& col_idx, const std::vector<double>& values, int m, int n);
        COO(const BaseMatrix& A);
        ~COO();
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
        std::vector<int> row_idx;
        std::vector<int> col_idx;
        std::vector<double> values;
        int m; // number of rows
        int n; // number of columns
};  

#endif
