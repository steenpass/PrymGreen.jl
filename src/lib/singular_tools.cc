#include "singular_tools.h"

/*
 * modified copy from Singular.jl (commit 0b7d433) to make use_cache and
 * use_tensor_trick available.
 */
std::tuple<void *, int> id_fres(ideal I, int n, std::string method,
        bool use_cache, bool use_tensor_trick, ring R)
{
    const ring origin = currRing;
    rChangeCurrRing(R);
    syStrategy s = syFrank(I, n, method.c_str(), use_cache, use_tensor_trick);
    rChangeCurrRing(origin);
    return std::make_tuple(reinterpret_cast<void *>(s->fullres), s->length);
}
