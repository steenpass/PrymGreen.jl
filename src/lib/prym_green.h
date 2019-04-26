#ifndef PRYM_GREEN_H
#define PRYM_GREEN_H

#include <stdlib.h>
#include <string.h>
#include <flint/ulong_extras.h>

#include "prym_green_types.h"
#include "binomial_coeffs.h"

msize_t multiply_matrix(arith_t **v_a, arith_t *A, arith_t* v_b, int g,
        arith_t c);
msize_t recurrence_sequence(arith_t **seq, arith_t *A, nvals_t n_values,
        arith_t* v, msize_t prym_green_size, msize_t index, int g, arith_t c);
msize_t kernel(arith_t **ker, arith_t *C, arith_t *A, nvals_t n_values,
        arith_t* v, msize_t prym_green_size, int g, arith_t c);
ulong mult_preinv_itest(ulong a, ulong b, ulong n, long N);

#endif   // PRYM_GREEN_H
