# 18-847A-Project

## Team Members
- Xiaobai Jiang (`xiaobaij`)
- Yuwei An (`yuweia`)
- Shuting Zhang (`shuting2`)

## Build Instructions

To run the tests, which include basic functionality tests for the matrix library and matrix decompositions, use:

```bash
make test
```
This will by default build the project in `release` mode, which is optimized for performance. If you want to build in `debug` mode, use:
```bash
make test BUILD=debug
```

To test SVD functionality specifically, use:
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