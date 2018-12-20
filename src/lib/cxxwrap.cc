#include "cxxwrap.h"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.method("return_msize_t", []() { return (msize_t)0; } );
    mod.method("return_nvals_t", []() { return (nvals_t)0; } );
    mod.method("return_entry_t", []() { return (entry_t)0; } );
    mod.method("return_arith_t", []() { return (arith_t)0; } );

    mod.method("id_fres", &id_fres);
    mod.method("rOrdStr", &rOrdStr_wrapper);
    mod.method("Singular_MaxBytesSystem", &Singular_MaxBytesSystem);
    mod.method("p_SubstPoly", &p_SubstPoly_wrapper);
}
