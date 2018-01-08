#include <flint/ulong_extras.h>

ulong preinv_test(ulong n)
{
    return n_preinvert_limb(n);
}

ulong mult_preinv_test(ulong a, ulong b, ulong n)
{
    ulong inv = n_preinvert_limb(n);
    return n_mulmod2_preinv(a, b, n, inv);
}

double mean(double a, double b)
{
    return (a+b) / 2;
}
