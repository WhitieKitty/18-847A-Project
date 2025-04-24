#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>    


class BaseMatrix {
    public:
        BaseMatrix();
        virtual ~BaseMatrix();
        virtual double get(int i, int j) const = 0;
        virtual void set(int i, int j, double value) = 0;
        virtual int getRows() const = 0;
        virtual int getCols() const = 0;
        virtual BaseMatrix* multiply(const BaseMatrix& other) const = 0;
        virtual BaseMatrix* add(const BaseMatrix& other) const = 0;
        virtual BaseMatrix* subtract(const BaseMatrix& other) const = 0;
        virtual BaseMatrix* transpose() const = 0;
        virtual void print() const = 0;
        virtual void show_matrix() const = 0;
};



#endif
