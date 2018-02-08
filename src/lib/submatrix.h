#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "kernel/ideals.h"

#pragma GCC diagnostic pop

#include "prym_green_types.h"

nvals_t check_matrix(entry_t **values_ptr, resolvente res, int g, msize_t size,
        msize_t limit, entry_t c, ring R);

#endif   // SUBMATRIX_H
