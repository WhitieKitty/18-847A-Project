# 18-847A-Project

## Team Members
- Xiaobai Jiang (`xiaobaij`)
- Yuwei An (`yuweia`)
- Shuting Zhang (`shuting2`)

## Build Instructions

To run the unit tests, which include basic functionality tests for the matrix library and matrix decompositions, use:

```bash
make test
```
To test SVD functionality:

```bash
make test_svd
```

## Dependencies
This project requires the following dependencies:
- C++11
- BLAS/LAPACK for Linux or WSL (`libblas` and `liblapacke`)
- Doxygen (for documentation generation)

## Example Usage

Example input/output tests are included in the `test` and `test_svd` targets. Their outpus demonstrate:
- Sparse matrix input and storage
- Eigenvalue and eigenvector computation
- SVD computation
- Matrix decomposition (QR, LU, Cholesky)

See `report.pdf` for a detailed explanation of the results.