#include "submatrix.h"

static int** init_binomial_coeffs(int g)
{
    int size = g-4;
    int **B = (int **)malloc(size * sizeof(int *));
    B[0] = (int *)calloc(size * size, sizeof(int));
    for (int i = 0; i < size; i++) {
        B[i] = (*B + size * i);
    }
    for (int i = 0; i < size; i++) {
        B[i][0] = 1;
        for (int j = 1; j <= i; j++) {
            B[i][j] = B[i-1][j-1]+B[i-1][j];
        }
    }
    return B;
}

static int binom(int n, int k, int **B)
{
    if (k < 0 || k > n) {
        return 0;
    }
    return B[n][k];
}

static long count_values_block(int h, int v, int f, int **B)
{
    if (h < 1 || v < 1 || f < 1) {
        return 0;
    }
    return binom(h+v+f-3, f-1, B);
}

static long count_values_row(int *hblocks, int n_hblocks, int g, int f,
	int **B)
{
    long n_values = 0;
    int v = g/2-f-1;
    for (int i = 0; i < n_hblocks; i++) {
	n_values += hblocks[i]*count_values_block(i+1, v, f, B);
	if (i != n_hblocks-1) {
	    n_values += hblocks[i]*count_values_block(i+1, v-1, f+1, B);
	}
    }
    return n_values;
}

static long count_values(int *hblocks, int n_hblocks, int g, int **B)
{
    long n_values = 0;
    n_values += 3*count_values_row(hblocks, n_hblocks, g, 2, B);
    n_values += g*count_values_row(hblocks, n_hblocks, g, 2, B);
    n_values += g*count_values_row(hblocks, n_hblocks, g, 3, B);
    return n_values;
}

int check_matrix(resolvente res, int g)
{
    int index = g/2-2;
//     module M = res[index];
    int **B = init_binomial_coeffs(g);

    /* define horizontal blocks */
    int n_hblocks = g/2-1;
    int hblocks[n_hblocks];
    for (int i = 0; i < n_hblocks-1; i++) {
        hblocks[i] = g/2+i-2;
    }
    hblocks[n_hblocks-1] = g-7;

    long n_values = count_values(hblocks, n_hblocks, g, B);
printf("n_values: %ld\n", n_values);

    return g;
}
