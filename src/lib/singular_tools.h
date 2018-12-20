#ifndef SINGULAR_TOOLS_H
#define SINGULAR_TOOLS_H

#include <tuple>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

#include "kernel/ideals.h"
#include "kernel/GBEngine/syz.h"

#pragma GCC diagnostic pop

std::tuple<void *, int> id_fres(ideal I, int n, std::string method,
        bool use_cache, bool use_tensor_trick, ring R);

#endif   // SINGULAR_TOOLS_H
