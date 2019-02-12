#ifndef SINGULAR_TOOLS_H
#define SINGULAR_TOOLS_H

#include <string>
#include <tuple>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "coeffs/longrat.h"
#include "kernel/ideals.h"
#include "kernel/GBEngine/syz.h"
#include "kernel/maps/gen_maps.h"

#pragma GCC diagnostic pop

std::tuple<void *, int> id_fres(ideal I, int n, std::string method,
        bool use_cache, bool use_tensor_trick, ring R);
std::string rOrdStr_wrapper(ring r);
long Singular_MaxBytesSystem();

#endif   // SINGULAR_TOOLS_H
