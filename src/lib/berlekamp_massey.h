#ifndef BERLEKAMP_MASSEY_H
#define BERLEKAMP_MASSEY_H

#include <string.h>
#include <flint/ulong_extras.h>

#include "prym_green_types.h"

msize_t berlekamp_massey(arith_t **lfsr, arith_t *seq, msize_t length,
        arith_t c);

#endif   // BERLEKAMP_MASSEY_H
