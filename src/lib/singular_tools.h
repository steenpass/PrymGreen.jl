#ifndef SINGULAR_TOOLS_H
#define SINGULAR_TOOLS_H

#include <tuple>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "kernel/ideals.h"
#include "kernel/GBEngine/syz.h"
#include "kernel/maps/gen_maps.h"

#pragma GCC diagnostic pop

std::tuple<void *, int> id_fres(ideal I, int n, std::string method,
        bool use_cache, bool use_tensor_trick, ring R);
std::string rOrdStr_wrapper(ring r);
long Singular_MaxBytesSystem();
poly p_SubstPoly_wrapper(poly p, int var, poly image, const ring preimage_r,
        const ring image_r);

#endif   // SINGULAR_TOOLS_H
