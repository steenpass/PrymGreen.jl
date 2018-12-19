#ifndef CXXWRAP_H
#define CXXWRAP_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include "jlcxx/jlcxx.hpp"

#pragma GCC diagnostic pop

#include "prym_green_types.h"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod);

#endif   // CXXWRAP_H
