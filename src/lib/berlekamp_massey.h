#ifndef BERLEKAMP_MASSEY_H
#define BERLEKAMP_MASSEY_H

#include "prym_green_types.h"

msize_t berlekamp_massey(arith_t **lfsr, arith_t *seq, msize_t length);

#endif   // BERLEKAMP_MASSEY_H
