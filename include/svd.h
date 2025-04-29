#ifndef SVD_H
#define SVD_H

#include "base_matrix.h"
#include <vector>

class SVD {
    public:
        static double lanczos_bidiag(const BaseMatrix& A, int k);
};


#endif