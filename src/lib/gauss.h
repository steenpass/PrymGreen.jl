#ifndef GAUSS_H
#define GAUSS_H

#include <string.h>
#include <flint/ulong_extras.h>

#include "prym_green_types.h"

msize_t gauss(entry_t *A, msize_t n_rows, msize_t n_cols, entry_t c);

#endif   // GAUSS_H
