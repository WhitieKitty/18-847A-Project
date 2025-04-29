#include "coo.h"

COO::COO(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
}

COO::~COO() {
}

void COO::print() const {
    std::cout << "Row Indices: ";
    for (int i = 0; i < row_idx.size(); i++) {
        std::cout << row_idx[i] << " "; 
    }
    std::cout << std::endl;
    std::cout << "Column Indices: ";
    for (int i = 0; i < col_idx.size(); i++) {
        std::cout << col_idx[i] << " ";
    }   
    std::cout << std::endl;
    std::cout << "Values: ";
    for (int i = 0; i < values.size(); i++) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl; 
}

void COO::show_matrix() const {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            bool found = false;
            for (int k = 0; k < values.size(); k++) {
                if (row_idx[k] == i && col_idx[k] == j) {
                    std::cout << values[k] << " ";
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
}

COO::COO(const std::vector<int>& row_idx, const std::vector<int>& col_idx, const std::vector<double>& values, int m, int n) {
    this->row_idx = row_idx;
    this->col_idx = col_idx;
    this->values = values;
    this->m = m;
    this->n = n;
}

BaseMatrix* COO::multiply(const BaseMatrix& other) const {
    const COO* result = dynamic_cast<const COO*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a COO matrix");
    }
    if (n != result->getRows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    std::vector<int> new_row_idx;
    std::vector<int> new_col_idx; 
    std::vector<double> new_values;

    for (size_t i = 0; i < values.size(); i++) {
        int row_a = row_idx[i];
        int col_a = col_idx[i];
        double val_a = values[i];

        for (size_t j = 0; j < result->values.size(); j++) {
            int row_b = result->row_idx[j];
            int col_b = result->col_idx[j];
            double val_b = result->values[j];

            if (col_a == row_b) {
                double product = val_a * val_b;
                
                bool found = false;
                for (size_t k = 0; k < new_values.size(); k++) {
                    if (new_row_idx[k] == row_a && new_col_idx[k] == col_b) {
                        new_values[k] += product;
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    new_row_idx.push_back(row_a);
                    new_col_idx.push_back(col_b);
                    new_values.push_back(product);
                }
            }
        }
    }

    return new COO(new_row_idx, new_col_idx, new_values, m, result->getCols());
}

BaseMatrix* COO::add(const BaseMatrix& other) const {
    const COO* result = dynamic_cast<const COO*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a COO matrix");
    }
    if (m != result->getRows() || n != result->getCols()) {
        throw std::invalid_argument("Matrix dimensions do not match for addition");
    }

    std::vector<int> new_row_idx;
    std::vector<int> new_col_idx;
    std::vector<double> new_values;

    size_t i = 0, j = 0;
    while (i < values.size() || j < result->values.size()) {
        if (i >= values.size()) {
            new_row_idx.push_back(result->row_idx[j]);
            new_col_idx.push_back(result->col_idx[j]);
            new_values.push_back(result->values[j]);
            j++;
        }
        else if (j >= result->values.size()) {
            new_row_idx.push_back(row_idx[i]);
            new_col_idx.push_back(col_idx[i]);
            new_values.push_back(values[i]);
            i++;
        }
        else {
            int row_a = row_idx[i];
            int col_a = col_idx[i];
            int row_b = result->row_idx[j];
            int col_b = result->col_idx[j];

            if (row_a < row_b || (row_a == row_b && col_a < col_b)) {
                new_row_idx.push_back(row_a);
                new_col_idx.push_back(col_a);
                new_values.push_back(values[i]);
                i++;
            }
            else if (row_b < row_a || (row_b == row_a && col_b < col_a)) {
                new_row_idx.push_back(row_b);
                new_col_idx.push_back(col_b);
                new_values.push_back(result->values[j]);
                j++;
            }
            else {
                new_row_idx.push_back(row_a);
                new_col_idx.push_back(col_a);
                new_values.push_back(values[i] + result->values[j]);
                i++;
                j++;
            }
        }
    }

    return new COO(new_row_idx, new_col_idx, new_values, m, n);
}

BaseMatrix* COO::subtract(const BaseMatrix& other) const {
    const COO* result = dynamic_cast<const COO*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a COO matrix");
    }
    if (m != result->m || n != result->n) {
        throw std::invalid_argument("Matrix dimensions do not match");
    }

    std::vector<int> new_row_idx;
    std::vector<int> new_col_idx;
    std::vector<double> new_values;

    size_t i = 0, j = 0;
    while (i < values.size() || j < result->values.size()) {
        if (i >= values.size()) {
            new_row_idx.push_back(result->row_idx[j]);
            new_col_idx.push_back(result->col_idx[j]);
            new_values.push_back(-result->values[j]);
            j++;
        }
        else if (j >= result->values.size()) {
            new_row_idx.push_back(row_idx[i]);
            new_col_idx.push_back(col_idx[i]);
            new_values.push_back(values[i]);
            i++;
        }
        else {
            int row_a = row_idx[i];
            int col_a = col_idx[i];
            int row_b = result->row_idx[j];
            int col_b = result->col_idx[j];

            if (row_a < row_b || (row_a == row_b && col_a < col_b)) {
                new_row_idx.push_back(row_a);
                new_col_idx.push_back(col_a);
                new_values.push_back(values[i]);
                i++;
            }
            else if (row_b < row_a || (row_b == row_a && col_b < col_a)) {
                new_row_idx.push_back(row_b);
                new_col_idx.push_back(col_b);
                new_values.push_back(-result->values[j]);
                j++;
            }
            else {
                new_row_idx.push_back(row_a);
                new_col_idx.push_back(col_a);
                new_values.push_back(values[i] - result->values[j]);
                i++;
                j++;
            }
        }
    }

    return new COO(new_row_idx, new_col_idx, new_values, m, n);
}

BaseMatrix* COO::transpose() const {
    std::vector<int> new_row_idx(col_idx);
    std::vector<int> new_col_idx(row_idx);
    std::vector<double> new_values(values);
    return new COO(new_row_idx, new_col_idx, new_values, n, m);
    
}

int COO::getRows() const {
    return m;
}

int COO::getCols() const {
    return n;
}

std::vector<double> COO::multiply(const std::vector<double>& other) const {
    if ((int)other.size() != n) {
        throw std::invalid_argument("Vector dimension does not match matrix columns");
    }

    std::vector<double> result(m, 0.0);
    for (size_t i = 0; i < values.size(); i++) {
        result[row_idx[i]] += values[i] * other[col_idx[i]];
    }
    return result;
}

void COO::set(int i, int j, double value) {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (size_t idx = 0; idx < row_idx.size(); idx++) {
        if (row_idx[idx] == i && col_idx[idx] == j) {
            values[idx] = value;
            return;
        }
    }

    row_idx.push_back(i);
    col_idx.push_back(j);
    values.push_back(value);
}

double COO::get(int i, int j) const {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (size_t idx = 0; idx < row_idx.size(); idx++) {
        if (row_idx[idx] == i && col_idx[idx] == j) {
            return values[idx];
        }
    }

    return 0.0;
}
