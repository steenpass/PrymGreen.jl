#include "prym_green.h"

msize_t recurrence_sequence(ulong **seq, entry_t *A, nvals_t n_values)
{
    int N = 10;
    *seq = (ulong *)malloc(N*sizeof(ulong));
    for (int i = 0; i < N; i++) {
        (*seq)[i] = 2*i;
    }
    return (msize_t)N;
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
