#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

/**
 * @class BaseMatrix
 * @brief Abstract base class for matrix operations.
 * 
 * This class defines the interface for a general-purpose matrix object.
 * Derived classes should implement the virtual methods for various matrix operations.
 */
class BaseMatrix {
public:
    /**
     * @brief Constructor.
     */
    BaseMatrix();

    /**
     * @brief Virtual destructor.
     */
    virtual ~BaseMatrix();

    /**
     * @brief Get the value at position (i, j).
     * @param i Row index.
     * @param j Column index.
     * @return The value at position (i, j).
     */
    virtual double get(int i, int j) const = 0;

    /**
     * @brief Set the value at position (i, j).
     * @param i Row index.
     * @param j Column index.
     * @param value The value to set.
     */
    virtual void set(int i, int j, double value) = 0;

    /**
     * @brief Get the number of rows.
     * @return The number of rows.
     */
    virtual int getRows() const = 0;

    /**
     * @brief Get the number of columns.
     * @return The number of columns.
     */
    virtual int getCols() const = 0;

    /**
     * @brief Multiply this matrix with another matrix.
     * @param other The other matrix.
     * @return A new matrix representing the product.
     */
    virtual BaseMatrix* multiply(const BaseMatrix& other) const = 0;

    /**
     * @brief Multiply this matrix with a vector.
     * @param other The vector.
     * @return A new vector representing the result.
     */
    virtual std::vector<double> multiply(const std::vector<double>& other) const = 0;

    /**
     * @brief Add another matrix to this matrix.
     * @param other The other matrix.
     * @return A new matrix representing the sum.
     */
    virtual BaseMatrix* add(const BaseMatrix& other) const = 0;

    /**
     * @brief Subtract another matrix from this matrix.
     * @param other The other matrix.
     * @return A new matrix representing the difference.
     */
    virtual BaseMatrix* subtract(const BaseMatrix& other) const = 0;

    /**
     * @brief Transpose this matrix.
     * @return A new matrix representing the transpose.
     */
    virtual BaseMatrix* transpose() const = 0;

    /**
     * @brief Print the matrix to standard output according to the format.
     */
    virtual void print() const = 0;

    /**
     * @brief Display the matrix content in a user-friendly format.
     */
    virtual void show_matrix() const = 0;
};

#endif // BASE_MATRIX_H
