# Sparse Matrix Library

A C++ library for handling sparse matrices with various storage formats and operations, including matrix decomposition and eigenvalue computation.

## Features

- Multiple sparse matrix storage formats:
  - COO (Coordinate Format)
  - CSR (Compressed Sparse Row)
  - CSC (Compressed Sparse Column)
- Matrix operations and decompositions:
  - Singular Value Decomposition (SVD)
  - Eigenvalue computation
  - Matrix decomposition utilities
- Matrix generation utilities
- Comprehensive test suite

## Prerequisites

- C++11 compatible compiler (clang++ recommended)
- LAPACK and BLAS libraries
  - On macOS: Accelerate framework
  - On Linux: liblapacke, liblapack, and libblas

## Building the Project

The project uses a Makefile-based build system. You can build the project in either debug or release mode.

### Debug Build
```bash
make BUILD=debug
```

### Release Build
```bash
make
```

### Running Tests
```bash
make test
make test_svd
```

### Documentation
To generate documentation using Doxygen:
```bash
make doc
```

## Project Structure

- `include/` - Header files
- `src/` - Source files
- `bin/` - Compiled binaries
- `obj/` - Object files
- `docs/` - Generated documentation

## Key Components

- `base_matrix.h/cpp` - Base class for sparse matrices
- `coo.h/cpp` - Coordinate format implementation
- `csr.h/cpp` - Compressed Sparse Row implementation
- `csc.h/cpp` - Compressed Sparse Column implementation
- `svd.h/cpp` - Singular Value Decomposition implementation
- `EigenSolver.h/cpp` - Eigenvalue computation
- `Decomposition.h/cpp` - Matrix decomposition utilities
- `matrix_generator.h/cpp` - Matrix generation utilities

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]