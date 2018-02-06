#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <cstdint>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "kernel/ideals.h"

#pragma GCC diagnostic pop

#define msize_t uint32_t   // for size of matrix and for module components
#define nvals_t uint32_t   // for number of values
#define entry_t uint16_t   // for entries and characteristic
#define null_entry ((entry_t)65535)   // represents not yet assigned entry

nvals_t check_matrix(entry_t **values_ptr, resolvente res, int g, msize_t size,
        msize_t limit, entry_t c, ring R);

#endif   // SUBMATRIX_H
