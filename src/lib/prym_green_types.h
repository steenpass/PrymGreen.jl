#ifndef PRYM_GREEN_TYPES_H
#define PRYM_GREEN_TYPES_H

#include <stdint.h>
#include <flint/flint.h>

#if FLINT_BITS != 64
#error "not implemented for FLINT_BITS != 64"
#endif

#define msize_t uint32_t   // for size of matrix and for module components
#define nvals_t uint32_t   // for number of values
#define entry_t uint32_t   // for entries and characteristic
#define null_entry ((entry_t)4294967295)   // represents not yet assigned entry
// uint64_t is ulong
#define arith_t uint64_t   // for modular arithmetic
#define ARITH_T_SUFFIX "ULL"   // suffix for literals of type arith_t

#endif   // PRYM_GREEN_TYPES_H
