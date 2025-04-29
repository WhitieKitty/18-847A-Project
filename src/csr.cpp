#include "csr.h"

CSR::CSR(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
}

CSR::~CSR() {
}

void CSR::print() const {
    std::cout << "Row Pointers: ";
    for (int i = 0; i < row_ptr.size(); i++) {
        std::cout << row_ptr[i] << " ";
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

void CSR::show_matrix() const {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            bool found = false;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                if (col_idx[k] == j) {
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

CSR::CSR(const std::vector<int>& row_ptr, const std::vector<int>& col_idx, const std::vector<double>& values, int m, int n) {
    this->row_ptr = row_ptr;
    this->col_idx = col_idx;
    this->values = values;
    this->m = m;
    this->n = n;
}

BaseMatrix* CSR::multiply(const BaseMatrix& other) const {
    const CSR* result = dynamic_cast<const CSR*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSR matrix");
    }
    if (n != result->getRows()) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    std::vector<int> new_row_ptr;
    std::vector<int> new_col_idx;
    std::vector<double> new_values;

    for (int i = 0; i < m; i++) {
        new_row_ptr.push_back(new_values.size());
        for (int j = 0; j < result->n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += get(i, k) * result->get(k, j);
            }
            if (sum != 0.0) {
                new_col_idx.push_back(j);
                new_values.push_back(sum);
            }
        }
    }
    new_row_ptr.push_back(new_values.size());

    return new CSR(new_row_ptr, new_col_idx, new_values, m, result->n);
}

std::vector<double> CSR::multiply(const std::vector<double>& other) const {
    if ((int)other.size() != n) {
        throw std::invalid_argument("Vector dimension does not match matrix columns");
    }

    std::vector<double> result(m, 0.0);
    for (int i = 0; i < m; i++) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
            result[i] += values[k] * other[col_idx[k]];
        }
    }
    return result;
}


BaseMatrix* CSR::add(const BaseMatrix& other) const {
    const CSR* result = dynamic_cast<const CSR*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSR matrix");
    }
    if (m != result->m || n != result->n) {
        throw std::invalid_argument("Matrix dimensions do not match for addition");
    }

    std::vector<int> new_row_ptr;
    std::vector<int> new_col_idx;
    std::vector<double> new_values;

    for (int i = 0; i < m; i++) {
        new_row_ptr.push_back(new_values.size());
        for (int j = 0; j < n; j++) {
            double sum = get(i, j) + result->get(i, j);
            if (sum != 0.0) {
                new_col_idx.push_back(j);
                new_values.push_back(sum);
            }
        }
    }
    new_row_ptr.push_back(new_values.size());

    return new CSR(new_row_ptr, new_col_idx, new_values, m, n);
}

BaseMatrix* CSR::subtract(const BaseMatrix& other) const {
    const CSR* result = dynamic_cast<const CSR*>(&other);
    if (result == nullptr) {
        throw std::invalid_argument("B matrix is not a CSR matrix");
    }
    if (m != result->m || n != result->n) {
        throw std::invalid_argument("Matrix dimensions do not match for subtraction");
    }

    std::vector<int> new_row_ptr;
    std::vector<int> new_col_idx;
    std::vector<double> new_values;

    for (int i = 0; i < m; i++) {
        new_row_ptr.push_back(new_values.size());
        for (int j = 0; j < n; j++) {
            double diff = get(i, j) - result->get(i, j);
            if (diff != 0.0) {
                new_col_idx.push_back(j);
                new_values.push_back(diff);
            }
        }
    }
    new_row_ptr.push_back(new_values.size());

    return new CSR(new_row_ptr, new_col_idx, new_values, m, n);
}

BaseMatrix* CSR::transpose() const {
    std::vector<int> new_row_ptr(n + 1, 0);
    std::vector<int> new_col_idx(values.size());
    std::vector<double> new_values(values.size());

    for (int i = 0; i < values.size(); i++) {
        new_row_ptr[col_idx[i] + 1]++;
    }

    for (int i = 1; i <= n; i++) {
        new_row_ptr[i] += new_row_ptr[i - 1];
    }

    std::vector<int> col_count(n, 0);
    for (int i = 0; i < m; i++) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            int col = col_idx[j];
            int pos = new_row_ptr[col] + col_count[col];
            new_col_idx[pos] = i;
            new_values[pos] = values[j];
            col_count[col]++;
        }
    }

    return new CSR(new_row_ptr, new_col_idx, new_values, n, m);
}

int CSR::getRows() const {
    return m;
}

int CSR::getCols() const {
    return n;
}

void CSR::set(int i, int j, double value) {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
        if (col_idx[k] == j) {
            values[k] = value;
            return;
        }
    }

    int pos = row_ptr[i + 1];
    values.insert(values.begin() + pos, value);
    col_idx.insert(col_idx.begin() + pos, j);
    for (int k = i + 1; k <= m; k++) {
        row_ptr[k]++;
    }
}

double CSR::get(int i, int j) const {
    if (i < 0 || i >= m || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of bounds");
    }

    for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
        if (col_idx[k] == j) {
            return values[k];
        }
    }

    return 0.0;
}