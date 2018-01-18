#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "polys/simpleideals.h"
#include "kernel/ideals.h"

#pragma GCC diagnostic pop

long check_matrix(int **values, resolvente res, int g, int64_t size,
        int limit);

#endif   // SUBMATRIX_H
