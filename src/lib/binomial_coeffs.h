#ifndef BINOMIAL_COEFFS_H
#define BINOMIAL_COEFFS_H

#include <stdlib.h>

#include "prym_green_types.h"

#define h_shift(h, v, f, B) binom(h+v+f-3, h-1, B);
#define v_shift(h, v, f, B) binom(h+v+f-3, v-1, B);

#ifdef __cplusplus
extern "C" {
#endif

int** init_binomial_coeffs(int g);
void clear_binomial_coeffs(int **B);
int binom(int n, int k, int **B);
msize_t prym_green_size(int g, int **B);
nvals_t count_values_block(int h, int v, int f, int **B);

#ifdef __cplusplus
}
#endif

#endif   // BINOMIAL_COEFFS_H
