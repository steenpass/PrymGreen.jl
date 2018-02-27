#ifndef PRYM_GREEN_H
#define PRYM_GREEN_H

#include <stdlib.h>
#include <string.h>
#include <flint/ulong_extras.h>

#include "prym_green_types.h"
#include "binomial_coeffs.h"

msize_t recurrence_sequence(arith_t **seq, entry_t *A, nvals_t n_values,
        entry_t* v, msize_t prym_green_size, msize_t index, int g, entry_t c);
ulong mult_preinv_itest(ulong a, ulong b, ulong n, long N);

#endif   // PRYM_GREEN_H
