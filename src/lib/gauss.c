#include "gauss.h"

msize_t gauss(entry_t *A_jl, msize_t n_rows, msize_t n_cols, entry_t c_jl)
{
    const arith_t c = (arith_t)c_jl;
    const arith_t c_inv = n_preinvert_limb(c);
    arith_t **A = (arith_t **)malloc(n_rows*sizeof(arith_t *));
    *A = (arith_t *)malloc(n_rows*n_cols*sizeof(arith_t));
    for (msize_t i = 1; i < n_rows; i++) {
        A[i] = A[i-1]+n_cols;
    }
    for (msize_t i = 0; i < n_rows; i++) {
        for (msize_t j = 0; j < n_cols; j++) {
            A[i][j] = (arith_t)A_jl[j*n_rows+i];
        }
    }
    msize_t rank = 0;
    msize_t current_row = 0;
    arith_t *tmp_row = (arith_t *)malloc(n_cols*sizeof(arith_t));
    for (msize_t j = 0; j < n_cols; j++) {
        msize_t i;
        for (i = current_row; i < n_rows; i++) {
            if (A[i][j] != 0) {
                break;
            }
        }
        if (i == n_rows) {
            continue;
        }
        rank++;
        if (i != current_row) {
            const size_t size = (n_cols-j)*sizeof(arith_t);
            memcpy(tmp_row, &A[i][j], size);
            memcpy(&A[i][j], &A[current_row][j], size);
            memcpy(&A[current_row][j], tmp_row, size);
        }
        if (A[current_row][j] != 1) {
            const arith_t inv = n_invmod(A[current_row][j], c);
            A[current_row][j] = 1;
            for (msize_t k = j+1; k < n_cols; k++) {
                A[current_row][k] = n_mulmod2_preinv(A[current_row][k], inv, c,
                        c_inv);
            }
        }
        for (msize_t i = current_row+1; i < n_rows; i++) {
            if (A[i][j] != 0) {
                for (msize_t k = j+1; k < n_cols; k++) {
                    const arith_t prod = n_mulmod2_preinv(A[i][j],
                            A[current_row][k], c, c_inv);
                    A[i][k] = n_submod(A[i][k], prod, c);
                }
                A[i][j] = 0;
            }
        }
        current_row++;
    }
    free(tmp_row);
    free(A[0]);
    free(A);
    return rank;
}
