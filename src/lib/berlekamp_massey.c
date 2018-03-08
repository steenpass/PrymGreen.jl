#include "berlekamp_massey.h"

static arith_t discrepancy(arith_t *C, arith_t *seq, msize_t N, msize_t L,
        arith_t c, arith_t inv)
{
    arith_t d = seq[N];
    for (msize_t i = 1; i <= L; i++) {
        arith_t prod = n_mulmod2_preinv(C[i], seq[N-i], c, inv);
        d = n_addmod(d, prod, c);
    }
    return d;
}

static void update_lfsr(arith_t *C, arith_t *B, msize_t L_B, arith_t d,
        arith_t b, msize_t x, arith_t c, arith_t inv)
{
    arith_t b_inv = n_invmod(b, c);
    arith_t prod = n_mulmod2_preinv(d, b_inv, c, inv);
    for (msize_t i = 0; i <= L_B; i++) {
        arith_t corr = n_mulmod2_preinv(prod, B[i], c, inv);
        C[x+i] = n_submod(C[x+i], corr, c);
    }
}

msize_t berlekamp_massey(arith_t **lfsr, arith_t *seq, msize_t length,
        arith_t c)
{
    // setup
    *lfsr = (arith_t *)calloc(length, sizeof(arith_t));
    arith_t inv = n_preinvert_limb(c);

    // step 1
    arith_t *C = *lfsr;
    C[0] = 1;
    msize_t L = 0;
    arith_t B[length];
    B[0] = 1;
    msize_t L_B = 0;   // length of B as LFSR
    arith_t b = 1;
    msize_t x = 1;

    // step 2
    arith_t d;
    arith_t T[length];
    for (msize_t N = 0; N < length; N++) {
        d = discrepancy(C, seq, N, L, c, inv);
        if (d == 0) {   // step 3
            x++;
        } else if (2*L > N) {   // step 4
            update_lfsr(C, B, L_B, d, b, x, c, inv);
            x++;
        } else {   // step 5
            memcpy(T, C, (L+1)*sizeof(arith_t));
            update_lfsr(C, B, L_B, d, b, x, c, inv);
            memcpy(B, T, (L+1)*sizeof(arith_t));
            L_B = L;
            L = N+1-L;
            b = d;
            x = 1;
        }
    }

    // clean-up and return
    *lfsr = (arith_t *)realloc(*lfsr, (L+1)*sizeof(arith_t));
    return L+1;
}
