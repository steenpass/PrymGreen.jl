#ifndef PRYM_GREEN_H
#define PRYM_GREEN_H

#include <flint/flint.h>
#include <flint/ulong_extras.h>

#include "prym_green_types.h"

#if FLINT_BITS != 64
#error "not implemented for FLINT_BITS != 64"
#endif

msize_t recurrence_sequence(ulong **seq, entry_t *A, nvals_t n_values,
        entry_t* v, msize_t prym_green_size, msize_t index, int g);
ulong mult_preinv_itest(ulong a, ulong b, ulong n, long N);

#endif   // PRYM_GREEN_H
