#include "prym_green.h"

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
