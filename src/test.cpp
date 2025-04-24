#include "matrix_generator.h"
#include "base_matrix.h"

int main() {
    MatrixGenerator mg;
    BaseMatrix* matrix = mg.generate_matrix("COO", 10, 10, 13);
    matrix->print();
    matrix->show_matrix();
    delete matrix;
}
