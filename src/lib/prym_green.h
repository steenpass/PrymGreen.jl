#ifndef PRYM_GREEN_H
#define PRYM_GREEN_H

#include <flint/flint.h>
#include <flint/ulong_extras.h>

#if FLINT_BITS != 64
#error "not implemented for FLINT_BITS != 64"
#endif

ulong mult_preinv_itest(ulong a, ulong b, ulong n, long N);

#endif   // PRYM_GREEN_H
