#include <flint/ulong_extras.h>
#include <flint/longlong.h>

#include "submatrix.h"

#if FLINT_BITS != 64
#error "not implemented for FLINT_BITS != 64"
#endif

int test_increase(int n)
{
    return increase(n);
}

ulong mult_preinv_test(ulong a, ulong b, ulong n, long N)
{
    ulong inv = n_preinvert_limb(n);
    ulong res;
    for (long i = 0; i <= N; i++) {
        res = n_mulmod2_preinv(a, b, n, inv);
        // res = n_addmod(a, b, n);
    }
    return res;
}
