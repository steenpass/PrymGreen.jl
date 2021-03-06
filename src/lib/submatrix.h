#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <stdlib.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "kernel/ideals.h"

#pragma GCC diagnostic pop

#include "prym_green_types.h"
#include "binomial_coeffs.h"

nvals_t check_matrix(entry_t **values_ptr, resolvente res, int g, msize_t size,
        msize_t limit, entry_t c, ring R);
msize_t dense_matrix(entry_t **M, resolvente res, int g, msize_t size,
        msize_t limit, ring R);

#endif   // SUBMATRIX_H
