#include "cxxwrap.h"

static nvals_t check_matrix_wrapper(void *values_ptr, void *res, int g,
        msize_t size, msize_t limit, entry_t c, ring R)
{
    return check_matrix(reinterpret_cast<entry_t **>(values_ptr),
            reinterpret_cast<resolvente>(res), g, size, limit, c, R);
}

static msize_t dense_matrix_wrapper(void *M, void *res, int g, msize_t size,
    msize_t limit, ring R)
{
    return dense_matrix(reinterpret_cast<entry_t **>(M),
            reinterpret_cast<resolvente>(res), g, size, limit, R);
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("return_msize_t", []() { return (msize_t)0; } );
    mod.method("return_nvals_t", []() { return (nvals_t)0; } );
    mod.method("return_entry_t", []() { return (entry_t)0; } );
    mod.method("return_arith_t", []() { return (arith_t)0; } );

    mod.method("id_fres", &id_fres);
    mod.method("rOrdStr", &rOrdStr_wrapper);
    mod.method("Singular_MaxBytesSystem", &Singular_MaxBytesSystem);

    mod.method("n_Int", [](number n, const coeffs r) { return n_Int(n, r); } );
    mod.method("p_LmCmp", &p_LmCmp);

    mod.method("check_matrix", &check_matrix_wrapper);
    mod.method("dense_matrix", &dense_matrix_wrapper);
}
