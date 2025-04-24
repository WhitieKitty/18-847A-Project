#include "csc.h"

CSC::CSC(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
}

CSC::~CSC() {
}

void CSC::print() const {
    std::cout << "Column Pointers: ";
    for (int i = 0; i < col_ptr.size(); i++) {
        std::cout << col_ptr[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Row Indices: ";
    for (int i = 0; i < row_idx.size(); i++) {
        std::cout << row_idx[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Values: ";
    for (int i = 0; i < values.size(); i++) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl;
}

void CSC::show_matrix() const {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            bool found = false;
            for (int k = col_ptr[j]; k < col_ptr[j + 1]; k++) {
                if (row_idx[k] == i) {
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

CSC::CSC(const std::vector<int>& col_ptr, const std::vector<int>& row_idx, const std::vector<double>& values, int m, int n) {
    this->col_ptr = col_ptr;
    this->row_idx = row_idx;
    this->values = values;
    this->m = m;
    this->n = n;
}

BaseMatrix* CSC::multiply(const BaseMatrix& other) const {
    const CSC* result = dynamic_cast<const CSC*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSC matrix");
    }
    if (n != result->getRows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    std::vector<int> new_col_ptr;
    std::vector<int> new_row_idx;
    std::vector<double> new_values;

    for (int j = 0; j < result->n; j++) {
        new_col_ptr.push_back(new_values.size());
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += get(i, k) * result->get(k, j);
            }
            if (sum != 0.0) {
                new_row_idx.push_back(i);
                new_values.push_back(sum);
            }
        }
    }
    new_col_ptr.push_back(new_values.size());

    return new CSC(new_col_ptr, new_row_idx, new_values, m, result->n);
}

BaseMatrix* CSC::add(const BaseMatrix& other) const {
    const CSC* result = dynamic_cast<const CSC*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSC matrix");
    }
    if (m != result->m || n != result->n) {
        throw std::invalid_argument("Matrix dimensions do not match for addition");
    }

    std::vector<int> new_col_ptr;
    std::vector<int> new_row_idx;
    std::vector<double> new_values;

    for (int j = 0; j < n; j++) {
        new_col_ptr.push_back(new_values.size());
        for (int i = 0; i < m; i++) {
            double sum = get(i, j) + result->get(i, j);
            if (sum != 0.0) {
                new_row_idx.push_back(i);
                new_values.push_back(sum);
            }
        }
    }
    new_col_ptr.push_back(new_values.size());

    return new CSC(new_col_ptr, new_row_idx, new_values, m, n);
}

BaseMatrix* CSC::subtract(const BaseMatrix& other) const {
    const CSC* result = dynamic_cast<const CSC*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSC matrix");
    }
    if (m != result->m || n != result->n) {
        throw std::invalid_argument("Matrix dimensions do not match for subtraction");
    }

    std::vector<int> new_col_ptr;
    std::vector<int> new_row_idx;
    std::vector<double> new_values;

    for (int j = 0; j < n; j++) {
        new_col_ptr.push_back(new_values.size());
        for (int i = 0; i < m; i++) {
            double diff = get(i, j) - result->get(i, j);
            if (diff != 0.0) {
                new_row_idx.push_back(i);
                new_values.push_back(diff);
            }
        }
    }
    new_col_ptr.push_back(new_values.size());

    return new CSC(new_col_ptr, new_row_idx, new_values, m, n);
}

BaseMatrix* CSC::transpose() const {
    std::vector<int> new_col_ptr(m + 1, 0);
    std::vector<int> new_row_idx(values.size());
    std::vector<double> new_values(values.size());

    for (int i = 0; i < values.size(); i++) {
        new_col_ptr[row_idx[i] + 1]++;
    }

    for (int i = 1; i <= m; i++) {
        new_col_ptr[i] += new_col_ptr[i - 1];
    }

    std::vector<int> row_count(m, 0);
    for (int j = 0; j < n; j++) {
        for (int i = col_ptr[j]; i < col_ptr[j + 1]; i++) {
            int row = row_idx[i];
            int pos = new_col_ptr[row] + row_count[row];
            new_row_idx[pos] = j;
            new_values[pos] = values[i];
            row_count[row]++;
        }
    }

    return new CSC(new_col_ptr, new_row_idx, new_values, n, m);
}

int CSC::getRows() const {
    return m;
}

int CSC::getCols() const {
    return n;
}

void CSC::set(int i, int j, double value) {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (int k = col_ptr[j]; k < col_ptr[j + 1]; k++) {
        if (row_idx[k] == i) {
            values[k] = value;
            return;
        }
    }

    int pos = col_ptr[j + 1];
    values.insert(values.begin() + pos, value);
    row_idx.insert(row_idx.begin() + pos, i);
    for (int k = j + 1; k <= n; k++) {
        col_ptr[k]++;
    }
}

double CSC::get(int i, int j) const {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (int k = col_ptr[j]; k < col_ptr[j + 1]; k++) {
        if (row_idx[k] == i) {
            return values[k];
        }
    }

    return 0.0;
}