#include "prym_green.h"

static void add_product(entry_t *v_a, entry_t *A, entry_t *v_b, int sign,
        entry_t c)
{
    ulong u_a = *v_a;
    ulong u_A = *A;
    ulong u_b = *v_b;
    ulong u_c = c;
    ulong inv = n_preinvert_limb(u_c);
    ulong res = n_mulmod2_preinv(u_A, u_b, u_c, inv);
    if (sign == -1) {
        res = u_c-res;
    }
    res = n_addmod(u_a, res, u_c);
    *v_a = (entry_t)res;
}

static void f_multiply(int h, int v, int sign, entry_t *v_a, entry_t *A,
        entry_t *v_b, entry_t c, int **B)
{
    int size = h_shift(h, v, 1, B);   // == v_shift(h, v, 1, B)
    for (int i = 0; i < size; i++) {
        add_product(&v_a[i], A, &v_b[i], sign, c);
    }
}

static void h_multiply(int v, int f, int sign, entry_t *v_a, entry_t *A,
        entry_t *v_b, entry_t c, int **B)
{
    int size = v_shift(1, v, f, B);
    for (int i = 0; i < size; i++) {
        add_product(&v_a[i], &A[size-i-1], v_b, sign, c);
        sign *= -1;
    }
}

static void v_multiply(int h, int f, int sign, entry_t *v_a, entry_t *A,
        entry_t *v_b, entry_t c, int **B)
{
    int size = h_shift(h, 1, f, B);
    for (int i = 0; i < size; i++) {
        add_product(v_a, &A[i], &v_b[i], sign, c);
    }
}

static void multiply_koszul_block(int h, int v, int f, int sign, entry_t *v_a,
        entry_t *A, entry_t *v_b, entry_t c, int **B)
{
    if (h < 1 || v < 1 || f < 1) {
        return;
    }
    if (f == 1) {
        return f_multiply(h, v, sign, v_a, A, v_b, c, B);
    } else if (h == 1 && f == 2) {
        return h_multiply(v, f, sign, v_a, A, v_b, c, B);
    } else if (v == 1) {
        return v_multiply(h, f, sign, v_a, A, v_b, c, B);
    }   // else:
    multiply_koszul_block(h-1, v, f, sign, v_a, A, v_b, c, B);
    v_b += h_shift(h-1, v, f, B);
    entry_t *f_A = &A[binom(h+v+f-4, f-1, B)];
    multiply_koszul_block(h, v, f-1, sign, v_a, f_A, v_b, c, B);
    v_a += v_shift(h, v, f-1, B);
    sign *= (f%2)*2-1;
    multiply_koszul_block(h, v-1, f, sign, v_a, A, v_b, c, B);
}

static void multiply_koszul_row(entry_t **v_a_iter, entry_t **A_iter,
        entry_t *v_b, int *hblocks, int n_hblocks, int g, int f, entry_t c,
        int **B)
{
    int v = g/2-f-1;
    entry_t *v_a = *v_a_iter;
    for (int i = 0; i < n_hblocks; i++) {
        for (int j = 0; j < hblocks[i]; j++) {
            multiply_koszul_block(i+1, v, f, 1, v_a, *A_iter, v_b, c, B);
            *A_iter += count_values_block(i+1, v, f, B);
            if (i != n_hblocks-1) {
                v_a += v_shift(i+1, v, f, B);
                multiply_koszul_block(i+1, v-1, f+1, 1, v_a, *A_iter, v_b, c,
                        B);
                *A_iter += count_values_block(i+1, v-1, f+1, B);
                v_a = *v_a_iter;
            }
            v_b += h_shift(i+1, v, f, B);
        }
    }
    *v_a_iter += v_shift(n_hblocks, v, f, B);
}

/*
 * let v_a = P * v_b where P is the Prym-Green matrix represented by A
 */
static void multiply_matrix(entry_t *v_a, entry_t *A, entry_t* v_b,
        int *hblocks, int n_hblocks, int g, entry_t c, int **B)
{
    entry_t **v_a_iter = (entry_t **)malloc(sizeof(entry_t *));
    *v_a_iter = v_a;
    entry_t **A_iter = (entry_t **)malloc(sizeof(entry_t *));
    *A_iter = A;
    for (int k = 0; k < 3; k++) {
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 2, c,
                B);
    }
    for (int k = 0; k < g; k++) {
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 2, c,
                B);
        multiply_koszul_row(v_a_iter, A_iter, v_b, hblocks, n_hblocks, g, 3, c,
                B);
    }
    free(v_a_iter);
    free(A_iter);
}

msize_t recurrence_sequence(ulong **seq, entry_t *A, nvals_t n_values,
        entry_t *v, msize_t prym_green_size, msize_t index, int g, entry_t c)
{
    int **B = init_binomial_coeffs(g);
    int **hblocks_ptr = (int **)malloc(sizeof(int *));
    int n_hblocks = init_horizontal_blocks(hblocks_ptr, g);
    size_t size_v = prym_green_size*sizeof(entry_t);
    entry_t *v_a = (entry_t *)malloc(size_v);
    entry_t *v_b = (entry_t *)malloc(size_v);
    memcpy(v_b, v, size_v);
    msize_t length = 2*prym_green_size;
    // this memory block will be handed over to the calling function:
    *seq = (ulong *)malloc(length*sizeof(ulong));
    (*seq)[0] = v[index];
    for (msize_t i = 1; i < length; i++) {
        memset(v_a, 0, size_v);
        multiply_matrix(v_a, A, v_b, *hblocks_ptr, n_hblocks, g, c, B);
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
