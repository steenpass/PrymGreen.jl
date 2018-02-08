#ifndef PRYM_GREEN_TYPES_H
#define PRYM_GREEN_TYPES_H

#include <stdint.h>

#define msize_t uint32_t   // for size of matrix and for module components
#define nvals_t uint32_t   // for number of values
#define entry_t uint16_t   // for entries and characteristic
#define null_entry ((entry_t)65535)   // represents not yet assigned entry

#endif   // PRYM_GREEN_TYPES_H
