#include "submatrix.h"

static nvals_t count_values_row(int *hblocks, int n_hblocks, int g, int f,
        int **B)
{
    nvals_t n_values = 0;
    int v = g/2-f-1;
    int i;
    for (i = 0; i < n_hblocks-1; i++) {
	n_values += (nvals_t)hblocks[i]*count_values_block(i+1, v, f, B);
	n_values += (nvals_t)hblocks[i]*count_values_block(i+1, v-1, f+1, B);
    }
    // i == n_hblocks-1
    n_values += (nvals_t)hblocks[i]*count_values_block(i+1, v, f, B);
    return n_values;
}

static nvals_t count_values(int *hblocks, int n_hblocks, int g, int **B)
{
    nvals_t n_values = 0;
    n_values += (nvals_t)3 * count_values_row(hblocks, n_hblocks, g, 2, B);
    n_values += (nvals_t)g * count_values_row(hblocks, n_hblocks, g, 2, B);
    n_values += (nvals_t)g * count_values_row(hblocks, n_hblocks, g, 3, B);
    return n_values;
}

static entry_t get_entry(poly *column, int sign, msize_t row, entry_t c)
{
    if (*column == NULL) {
        return (entry_t)0;
    }
    msize_t comp = (msize_t)pGetComp(*column);
    if (comp < row) {
        return null_entry;   // error: wrong entry
    }
    if (comp > row) {
        return (entry_t)0;
    }
    // else: comp == row
    entry_t m = (entry_t)(long)pGetCoeff(*column);
    pIter(*column);
    if (sign == -1) {   // assume m != 0
        m = c-m;
    }
    return m;
}

static bool check_entry(poly *column, entry_t *value, int sign, msize_t row,
        entry_t c)
{
    entry_t m = get_entry(column, sign, row, c);
    if (m == null_entry) {
        return false;
    }
    if (*value == null_entry) {
        *value = m;
        return true;
    }
    return (*value == m);
}

static bool f_check(int h, int v, int sign, poly *columns, entry_t *values,
        msize_t row, entry_t c, int **B)
{
    int size = h_shift(h, v, 1, B);   // == v_shift(h, v, 1, B)
    for (int i = 0; i < size; i++) {
        if (!check_entry(&columns[i], values, sign, row, c)) {
            return false;
        }
        row++;
    }
    return true;
}

static bool h_check(int v, int f, int sign, poly *columns, entry_t *values,
        msize_t row, entry_t c, int **B)
{
    int size = v_shift(1, v, f, B);
    for (int i = 0; i < size; i++) {
        if (!check_entry(columns, &values[size-i-1], sign, row, c)) {
            return false;
        }
        row++;
        sign *= -1;
    }
    return true;
}

static bool v_check(int h, int f, int sign, poly *columns, entry_t *values,
        msize_t row, entry_t c, int **B)
{
    int size = h_shift(h, 1, f, B);
    for (int i = 0; i < size; i++) {
        if (!check_entry(&columns[i], &values[i], sign, row, c)) {
            return false;
        }
    }
    return true;
}

static bool check_koszul_block(int h, int v, int f, int sign, poly *columns,
    entry_t *values, msize_t row, entry_t c, int **B)
{
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
        return false;
    }
    columns += h_shift(h-1, v, f, B);
    entry_t *f_values = &values[binom(h+v+f-4, f-1, B)];
    if (!check_koszul_block(h, v, f-1, sign, columns, f_values, row, c, B)) {
        return false;
    }
    row += (msize_t)v_shift(h, v, f-1, B);
    sign *= (f%2)*2-1;
    if (!check_koszul_block(h, v-1, f, sign, columns, values, row, c, B)) {
        return false;
    }
    return true;
}

static bool check_koszul_row(int *hblocks, int n_hblocks, int g, int f,
        poly *columns, entry_t **values_iter, msize_t *row_ptr, entry_t c,
        int **B)
{
    int v = g/2-f-1;
    msize_t row = *row_ptr;
    for (int i = 0; i < n_hblocks; i++) {
        for (int j = 0; j < hblocks[i]; j++) {
            if (!check_koszul_block(i+1, v, f, 1, columns, *values_iter, row,
                        c, B)) {
                return false;
            }
            *values_iter += count_values_block(i+1, v, f, B);
            if (i != n_hblocks-1) {
                row += (msize_t)v_shift(i+1, v, f, B);
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
    *row_ptr += (msize_t)v_shift(n_hblocks, v, f, B);
    return true;
}

static bool check_entries(entry_t **values_ptr, poly *columns, int *hblocks,
        int n_hblocks, int g, msize_t size, msize_t limit, entry_t c, int **B)
{
    /* check the rows of blocks */
    entry_t **values_iter = (entry_t **)malloc(sizeof(entry_t *));
    *values_iter = *values_ptr;
    msize_t row = limit+1;
    // three rows of blocks with f = 2
    for (int k = 0; k < 3; k++) {
        if (!check_koszul_row(hblocks, n_hblocks, g, 2, columns, values_iter,
                    &row, c, B)) {
            free(values_iter);
            return false;
        }
    }
    // g double rows of blocks with f = 2 and f = 3
    for (int k = 0; k < g; k++) {
        if (!check_koszul_row(hblocks, n_hblocks, g, 2, columns, values_iter,
                    &row, c, B)) {
            free(values_iter);
            return false;
        }
        if (!check_koszul_row(hblocks, n_hblocks, g, 3, columns, values_iter,
                    &row, c, B)) {
            free(values_iter);
            return false;
        }
    }
    free(values_iter);

    /* check that there are no entries left in any of the columns */
    for (msize_t i = 0; i < size; i++) {
        if (columns[i] != NULL) {
            return false;
        }
    }
    return true;
}

static nvals_t check_matrix_currRing(entry_t **values_ptr, resolvente res,
        int g, msize_t size, msize_t limit, entry_t c, int **B)
{
    /* check size */
    if (size != prym_green_size(g, B)) {
        fprintf(stderr, "error: Prym-Green size does not fit\n");
        return 0;   // error
    }

    /* define horizontal blocks */
    int **hblocks_ptr = (int **)malloc(sizeof(int *));
    int n_hblocks = init_horizontal_blocks(hblocks_ptr, g);

    /* define array of columns */
    poly *columns = (poly *)malloc(size*sizeof(poly));
    int index = g/2-3;
    for (msize_t i = 0; i < size; i++) {
        columns[i] = res[index]->m[i];
        while (columns[i] != NULL && (msize_t)pGetComp(columns[i]) <= limit) {
            pIter(columns[i]);
        }
    }

    /* define array of values */
    nvals_t n_values = count_values(*hblocks_ptr, n_hblocks, g, B);
    // this memory block will be handed over to the calling function:
    *values_ptr = (entry_t *)malloc(n_values*sizeof(entry_t));
    for (nvals_t i = 0; i < n_values; i++) {
        (*values_ptr)[i] = null_entry;
    }

    /* check entries and return */
    if (!check_entries(values_ptr, columns, *hblocks_ptr, n_hblocks, g, size,
            limit, c, B)) {
        fprintf(stderr, "error: matrix does not admit Prym-Green structure\n");
        n_values = 0;   // error
    }
    free(columns);
    clear_horizontal_blocks(hblocks_ptr);
    free(hblocks_ptr);
    return n_values;
}

nvals_t check_matrix(entry_t **values_ptr, resolvente res, int g, msize_t size,
        msize_t limit, entry_t c, ring R)
{
    const ring R_orig = currRing;
    rChangeCurrRing(R);
    int **B = init_binomial_coeffs(g);
    nvals_t n_values = check_matrix_currRing(values_ptr, res, g, size, limit,
            c, B);
    clear_binomial_coeffs(B);
    rChangeCurrRing(R_orig);
    return n_values;
}

/*
 * debug functions
 */

msize_t dense_matrix(ulong **M, resolvente res, int g, msize_t size,
        msize_t limit, ring R)
{
    const ring R_orig = currRing;
    rChangeCurrRing(R);
    poly *columns = (poly *)malloc(size*sizeof(poly));
    int index = g/2-3;
    for (msize_t i = 0; i < size; i++) {
        columns[i] = res[index]->m[i];
        while (columns[i] != NULL && (msize_t)pGetComp(columns[i]) <= limit) {
            pIter(columns[i]);
        }
    }
    *M = (ulong *)calloc(size*size, sizeof(ulong));
    for (msize_t i = 1; i < size; i++) {
        M[i] = M[i-1]+size;
    }
    for (msize_t i = 0; i < size; i++) {
        while (columns[i] != NULL) {
            msize_t comp = (msize_t)pGetComp(columns[i]);
            if (comp <= limit || comp > limit+size) {
                fprintf(stderr, "error: wrong component\n");
                free(columns);
                rChangeCurrRing(R_orig);
                return 0;   // error
            }
            // julia arrays are stored in column-major order
            M[i][comp-limit-1] = (ulong)(long)pGetCoeff(columns[i]);
            pIter(columns[i]);
        }
    }
    free(columns);
    rChangeCurrRing(R_orig);
    return size;
}
