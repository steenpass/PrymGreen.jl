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

static void clear_binomial_coeffs(int **B)
{
    free(B[0]);
    free(B);
}

static int binom(int n, int k, int **B)
{
    if (k < 0 || k > n) {
        return 0;
    }
    return B[n][k];
}

#define h_shift(h, v, f, B) binom(h+v+f-3, h-1, B)
#define v_shift(h, v, f, B) binom(h+v+f-3, v-1, B)

static int64_t prym_green_size(int *hblocks, int n_hblocks, int g, int **B)
{
    long size = 0;
    // sum up h_shift's of all blocks in the first row
    int f = 2;
    int v = g/2-f-1;
    for (int i = 0; i < n_hblocks; i++) {
	size += hblocks[i]*h_shift(i+1, v, f, B);
    }
    return size;
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

static bool check_entry(poly *column, int *value, int sign, int row, int c)
{
    if (*column == NULL || pGetComp(*column) != row) {
        return false;
    }
    int m = (int)(long)pGetCoeff(*column);
    pIter(*column);
    if (sign == -1 && m != 0) m = c-m;
// printf("checking entry %d ?= %d\n", m, *value);
    if (*value == -1) {
        *value = m;
        return true;
    } else {
        return (*value == m);
    }
}

static bool f_check(int h, int v, int sign, poly *columns, int *values,
        int row, int c, int **B)
{
// printf("f_check with (h, v, sign) = (%d, %d, %d)\n", h, v, sign);
    int size = h_shift(h, v, 1, B);   // == v_shift(h, v, 1, B)
    for (int i = 0; i < size; i++) {
        if (!check_entry(&columns[i], values, sign, row+i, c)) {
            return false;
        }
    }
    return true;
}

static bool h_check(int v, int f, int sign, poly *columns, int *values,
        int row, int c, int **B)
{
// printf("h_check with (v, f, sign) = (%d, %d, %d)\n", v, f, sign);
    int size = v_shift(1, v, f, B);
    for (int i = 0; i < size; i++) {
        if (!check_entry(columns, &values[size-i-1], sign, row+i, c)) {
            return false;
        }
        sign *= -1;
    }
    return true;
}

static bool v_check(int h, int f, int sign, poly *columns, int *values,
        int row, int c, int **B)
{
// printf("v_check with (h, f, sign) = (%d, %d, %d)\n", h, f, sign);
    int size = h_shift(h, 1, f, B);
    for (int i = 0; i < size; i++) {
        if (!check_entry(&columns[i], &values[i], sign, row, c)) {
            return false;
        }
    }
    return true;
}

static bool check_koszul_block(int h, int v, int f, int sign, poly *columns,
    int *values, int row, int c, int **B)
{
// printf("checking block with (h, v, f, sign) = (%d, %d, %d, %d)\n",
//          h, v, f, sign);
    if (h < 1 || v < 1 || f < 1) {
        return true;
    }
    if (f == 1) {
        return f_check(h, v, sign, columns, values, row, c, B);
    } else if (h == 1 && f == 2) {
        return h_check(v, f, sign, columns, values, row, c, B);
    } else if (v == 1) {
        return v_check(h, f, sign, columns, values, row, c, B);
    }   // else:
    if (!check_koszul_block(h-1, v, f, sign, columns, values, row, c, B)) {
        return 0;
    }
    columns += h_shift(h-1, v, f, B);
    if (!check_koszul_block(h, v, f-1, sign, columns,
                &values[binom(h+v+f-4, f-1, B)], row, c, B)) {
        return false;
    }
    row += v_shift(h, v, f-1, B);
    if (!check_koszul_block(h, v-1, f, ((f%2)*2-1)*sign, columns, values, row,
                c, B)) {
        return false;
    }
    return true;
}


static bool check_koszul_row(int hblocks[], int n_hblocks, int g, int f,
        poly *columns, int **values_iter, int *row_ptr, int c, int **B)
{
// printf("checking row with g = %d, f = %d\n", g, f);
    int v = g/2-f-1;
    int row = *row_ptr;
    for (int i = 0; i < n_hblocks; i++) {
        for (int j = 0; j < hblocks[i]; j++) {
            if (!check_koszul_block(i+1, v, f, 1, columns, *values_iter, row,
                        c, B)) {
                return false;
            }
            *values_iter += count_values_block(i+1, v, f, B);
            if (i != n_hblocks-1) {
                row += v_shift(i+1, v, f, B);
                if (!check_koszul_block(i+1, v-1, f+1, 1, columns,
                            *values_iter, row, c, B)) {
                    return false;
                }
                *values_iter += count_values_block(i+1, v-1, f+1, B);
                row = *row_ptr;
            }
            columns += h_shift(i+1, v, f, B);
        }
    }
    *row_ptr += v_shift(n_hblocks, v, f, B);
    return true;
}

static long check_matrix_currRing(int **values_ptr, resolvente res, int g,
        int64_t size, int limit, int c, int **B)
{
    int index = g/2-3;

    /* define horizontal blocks */
    int n_hblocks = g/2-1;
    int hblocks[n_hblocks];
    for (int i = 0; i < n_hblocks-1; i++) {
        hblocks[i] = g/2+i-2;
    }
    hblocks[n_hblocks-1] = g-7;

    if (size != prym_green_size(hblocks, n_hblocks, g, B)) {
        fprintf(stderr, "matrix not square, returning 0\n");
        return 0;   // error
    }
    poly *columns = (poly *)malloc(size*sizeof(poly));
    for (long i = 0; i < size; i++) {
        columns[i] = res[index]->m[i];
        while (columns[i] != NULL && pGetComp(columns[i]) <= limit) {
            pIter(columns[i]);
        }
    }
    long n_values = count_values(hblocks, n_hblocks, g, B);
    // this memory block will be handed over to the calling function:
    *values_ptr = (int *)malloc(n_values*sizeof(int));
    for (long i = 0; i < n_values; i++) {
        (*values_ptr)[i] = -1;
    }

    /* check entries */
    int **values_iter = (int **)malloc(sizeof(int *));
    *values_iter = *values_ptr;
    int row = limit+1;
    char errmsg[]
        = "matrix does not admit prym green structure, returning 0\n";
    for (int k = 0; k < 3; k++) {
        if (!check_koszul_row(hblocks, n_hblocks, g, 2, columns, values_iter,
                    &row, c, B)) {
            fprintf(stderr, errmsg);
            return 0;
        }
    }
    for (int k = 0; k < g; k++) {
        if (!check_koszul_row(hblocks, n_hblocks, g, 2, columns, values_iter,
                    &row, c, B)) {
            fprintf(stderr, errmsg);
            return 0;
        }
        if (!check_koszul_row(hblocks, n_hblocks, g, 3, columns, values_iter,
                    &row, c, B)) {
            fprintf(stderr, errmsg);
            return 0;
        }
    }
    for (long i = 0; i < size; i++) {
        if (columns[i] != NULL) {
            fprintf(stderr, errmsg);
            return 0;
        }
    }
    free(values_iter);
    free(columns);
    return n_values;
}

long check_matrix(int **values_ptr, resolvente res, int g, int64_t size,
        int limit, int c, ring R)
{
    const ring R_orig = currRing;
    rChangeCurrRing(R);
    int **B = init_binomial_coeffs(g);
    long n_values = check_matrix_currRing(values_ptr, res, g, size, limit, c,
            B);
    clear_binomial_coeffs(B);
    rChangeCurrRing(R_orig);
    return n_values;
}
