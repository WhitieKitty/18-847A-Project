#include "matrix_generator.h"

MatrixGenerator::MatrixGenerator() {
}

MatrixGenerator::~MatrixGenerator() {
}

BaseMatrix* MatrixGenerator::generate_matrix(const std::string& format, int m, int n, int density) {
    if (density < 0 || density > 100) {
        throw std::invalid_argument("Density must be between 0 and 100");
    }

    int total_elements = m * n;
    int non_zero_elements = (total_elements * density) / 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> row_dis(0, m - 1);
    std::uniform_int_distribution<> col_dis(0, n - 1);

    std::vector<int> row_indices;
    std::vector<int> col_indices;
    std::vector<double> values;

    for (int i = 0; i < non_zero_elements; i++) {
        int row = row_dis(gen);
        int col = col_dis(gen);
        double value = dis(gen);

        bool exists = false;
        for (size_t j = 0; j < row_indices.size(); j++) {
            if (row_indices[j] == row && col_indices[j] == col) {
                exists = true;
                break;
            }
        }

        if (!exists) {
            row_indices.push_back(row);
            col_indices.push_back(col);
            values.push_back(value);
        } else {
            i--;
        }
    }

    if (format == "CSR") {
        std::vector<int> row_ptr(m + 1, 0);
        std::vector<int> col_idx;
        std::vector<double> csr_values;

        for (size_t i = 0; i < row_indices.size(); i++) {
            row_ptr[row_indices[i] + 1]++;
        }

        for (int i = 1; i <= m; i++) {
            row_ptr[i] += row_ptr[i - 1];
        }

        std::vector<std::pair<int, int>> sorted_indices;
        for (size_t i = 0; i < row_indices.size(); i++) {
            sorted_indices.push_back({row_indices[i], col_indices[i]});
        }

        std::sort(sorted_indices.begin(), sorted_indices.end());

        col_idx.resize(row_indices.size());
        csr_values.resize(row_indices.size());

        for (size_t i = 0; i < sorted_indices.size(); i++) {
            col_idx[i] = sorted_indices[i].second;
            csr_values[i] = values[i];
        }

        return new CSR(row_ptr, col_idx, csr_values, m, n);
    } else if (format == "CSC") {
        std::vector<int> col_ptr(n + 1, 0);
        std::vector<int> row_idx;
        std::vector<double> csc_values;

        for (size_t i = 0; i < col_indices.size(); i++) {
            col_ptr[col_indices[i] + 1]++;
        }

        for (int i = 1; i <= n; i++) {
            col_ptr[i] += col_ptr[i - 1];
        }

        std::vector<std::pair<int, int>> sorted_indices;
        for (size_t i = 0; i < col_indices.size(); i++) {
            sorted_indices.push_back({col_indices[i], row_indices[i]});
        }

        std::sort(sorted_indices.begin(), sorted_indices.end());

        row_idx.resize(col_indices.size());
        csc_values.resize(col_indices.size());

        for (size_t i = 0; i < sorted_indices.size(); i++) {
            row_idx[i] = sorted_indices[i].second;
            csc_values[i] = values[i];
        }

        return new CSC(col_ptr, row_idx, csc_values, m, n);
    } else if (format == "COO") {
        return new COO(row_indices, col_indices, values, m, n);
    } else {
        throw std::invalid_argument("Invalid matrix format");
    }
} 

BaseMatrix* MatrixGenerator::generate_spd_matrix(int n) {
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_val(0.1, 1.0);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = distribution_val(generator);

    std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                ATA[i][j] += A[k][i] * A[k][j];

    std::vector<int> row_idx, col_idx;
    std::vector<double> values;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (ATA[i][j] != 0.0) {
                row_idx.push_back(i);
                col_idx.push_back(j);
                values.push_back(ATA[i][j]);
            }

    return new COO(row_idx, col_idx, values, n, n);
}
