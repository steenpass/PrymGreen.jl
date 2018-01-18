#include "polys/simpleideals.h"

int increase2(int n, ideal I)
{
    return ++n;
}

extern "C" int increase(int);

int increase(int n)
{
    ideal I;
    return increase2(n, I);
}
