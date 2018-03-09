#include "prym_green.h"

static void add_product(arith_t *v_a, arith_t *A, arith_t *v_b, int sign,
        arith_t c, arith_t inv)
{
    arith_t res = n_mulmod2_preinv(*A, *v_b, c, inv);
    if (sign == -1 && res != 0) {
        res = c-res;
    }
    *v_a = n_addmod(*v_a, res, c);
}

static void f_multiply(int h, int v, int sign, arith_t *v_a, arith_t *A,
        arith_t *v_b, arith_t c, arith_t inv, int **B)
{
    int size = h_shift(h, v, 1, B);   // == v_shift(h, v, 1, B)
    for (int i = 0; i < size; i++) {
        add_product(&v_a[i], A, &v_b[i], sign, c, inv);
    }
}

static void h_multiply(int v, int f, int sign, arith_t *v_a, arith_t *A,
        arith_t *v_b, arith_t c, arith_t inv, int **B)
{
    int size = v_shift(1, v, f, B);
    for (int i = 0; i < size; i++) {
        add_product(&v_a[i], &A[size-i-1], v_b, sign, c, inv);
        sign *= -1;
    }
}

static void v_multiply(int h, int f, int sign, arith_t *v_a, arith_t *A,
        arith_t *v_b, arith_t c, arith_t inv, int **B)
{
    int size = h_shift(h, 1, f, B);
    for (int i = 0; i < size; i++) {
        add_product(v_a, &A[i], &v_b[i], sign, c, inv);
    }
}

static void multiply_koszul_block(int h, int v, int f, int sign, arith_t *v_a,
        arith_t *A, arith_t *v_b, arith_t c, arith_t inv, int **B)
{
    if (h < 1 || v < 1 || f < 1) {
        return;
    }
    if (f == 1) {
        return f_multiply(h, v, sign, v_a, A, v_b, c, inv, B);
    } else if (h == 1 && f == 2) {
        return h_multiply(v, f, sign, v_a, A, v_b, c, inv, B);
    } else if (v == 1) {
        return v_multiply(h, f, sign, v_a, A, v_b, c, inv, B);
    }   // else:
    multiply_koszul_block(h-1, v, f, sign, v_a, A, v_b, c, inv, B);
    v_b += h_shift(h-1, v, f, B);
    arith_t *f_A = &A[binom(h+v+f-4, f-1, B)];
    multiply_koszul_block(h, v, f-1, sign, v_a, f_A, v_b, c, inv, B);
    v_a += v_shift(h, v, f-1, B);
    sign *= (f%2)*2-1;
    multiply_koszul_block(h, v-1, f, sign, v_a, A, v_b, c, inv, B);
}

static void multiply_koszul_row(arith_t **v_a_iter, arith_t **A_iter,
        arith_t *v_b, int *hblocks, int n_hblocks, int g, int f, arith_t c,
        arith_t inv, int **B)
{
    int v = g/2-f-1;
    arith_t *v_a = *v_a_iter;
    for (int i = 0; i < n_hblocks; i++) {
        for (int j = 0; j < hblocks[i]; j++) {
            multiply_koszul_block(i+1, v, f, 1, v_a, *A_iter, v_b, c, inv, B);
            *A_iter += count_values_block(i+1, v, f, B);
            if (i != n_hblocks-1) {
                v_a += v_shift(i+1, v, f, B);
                multiply_koszul_block(i+1, v-1, f+1, 1, v_a, *A_iter, v_b, c,
                        inv, B);
                *A_iter += count_values_block(i+1, v-1, f+1, B);
                v_a = *v_a_iter;
            }
            v_b += h_shift(i+1, v, f, B);
        }
    }
    *v_a_iter += v_shift(n_hblocks, v, f, B);
}

static void multiply_matrix_loop(arith_t *v_a, arith_t *A, arith_t* v_b,
        int *hblocks, int n_hblocks, int g, arith_t c, arith_t inv, int **B)
{
    arith_t **v_a_iter = (arith_t **)malloc(sizeof(arith_t *));
    *v_a_iter = v_a;
    arith_t **A_iter = (arith_t **)malloc(sizeof(arith_t *));
    *A_iter = A;
    for (int k = 0; k < 3; k++) {
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 2, c,
                inv, B);
    }
    for (int k = 0; k < g; k++) {
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 2, c,
                inv, B);
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 3, c,
                inv, B);
    }
    free(v_a_iter);
    free(A_iter);
}

/*
 * let v_a = P * v_b where P is the Prym-Green matrix represented by A
 */
msize_t multiply_matrix(arith_t **v_a, arith_t *A, arith_t* v_b, int g,
        arith_t c)
{
    int **B = init_binomial_coeffs(g);
    int **hblocks_ptr = (int **)malloc(sizeof(int *));
    int n_hblocks = init_horizontal_blocks(hblocks_ptr, g);
    arith_t inv = n_preinvert_limb(c);
    msize_t size = prym_green_size(g, B);
    *v_a = (arith_t *)calloc(size, sizeof(arith_t));
    multiply_matrix_loop(*v_a, A, v_b, *hblocks_ptr, n_hblocks, g, c, inv, B);
    clear_horizontal_blocks(hblocks_ptr);
    free(hblocks_ptr);
    clear_binomial_coeffs(B);
    return size;
}

msize_t recurrence_sequence(arith_t **seq, arith_t *A, nvals_t n_values,
        arith_t *v, msize_t prym_green_size, msize_t index, int g, arith_t c)
{
    int **B = init_binomial_coeffs(g);
    int **hblocks_ptr = (int **)malloc(sizeof(int *));
    int n_hblocks = init_horizontal_blocks(hblocks_ptr, g);
    arith_t inv = n_preinvert_limb(c);
    size_t size_v = prym_green_size*sizeof(arith_t);
    arith_t *v_a = (arith_t *)malloc(size_v);
    arith_t *v_b = (arith_t *)malloc(size_v);
    memcpy(v_b, v, size_v);
    msize_t length = 2*prym_green_size;
    // this memory block will be handed over to the calling function:
    *seq = (arith_t *)malloc(length*sizeof(arith_t));
    (*seq)[0] = v[index];
    for (msize_t i = 1; i < length; i++) {
        memset(v_a, 0, size_v);
        multiply_matrix_loop(v_a, A, v_b, *hblocks_ptr, n_hblocks, g, c, inv,
                B);
        (*seq)[i] = v_a[index];
        memcpy(v_b, v_a, size_v);
    }
    free(v_a);
    free(v_b);
    clear_horizontal_blocks(hblocks_ptr);
    free(hblocks_ptr);
    clear_binomial_coeffs(B);
    return length;
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
