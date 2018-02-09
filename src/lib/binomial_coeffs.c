#include "binomial_coeffs.h"

int** init_binomial_coeffs(int g)
{
    int size = g-3;
    int **B = (int **)malloc(size * sizeof(int *));
    *B = (int *)malloc((size*(size+1)/2) * sizeof(int));
    B[0][0] = 1;
    for (int n = 1; n < size; n++) {
        B[n] = B[n-1]+n;
        B[n][0] = 1;
        B[n][n] = 1;
        for (int k = 1; k < n; k++) {
            B[n][k] = B[n-1][k-1]+B[n-1][k];
        }
    }
    return B;
}

void clear_binomial_coeffs(int **B)
{
    free(B[0]);
    free(B);
}

int binom(int n, int k, int **B)
{
    if (k < 0 || k > n) {
        return 0;
    }
    return B[n][k];
}

msize_t prym_green_size(int g, int **B)
{
    // sum up v_shift's of all blocks in the last column
    msize_t size = (msize_t)3*(msize_t)binom(g-5, g/2-4, B);
    size += (msize_t)g*(msize_t)binom(g-4, g/2-4, B);
    return size;
}

nvals_t count_values_block(int h, int v, int f, int **B)
{
    if (h < 1 || v < 1 || f < 1) {
        return (nvals_t)0;
    }
    return (nvals_t)binom(h+v+f-3, f-1, B);
}
