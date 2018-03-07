#include "berlekamp_massey.h"

msize_t berlekamp_massey(arith_t **lfsr, arith_t *seq, msize_t length)
{
    *lfsr = (arith_t *)calloc(length, sizeof(arith_t));
    return length;
}
