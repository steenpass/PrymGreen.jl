#ifndef CXXWRAP_H
#define CXXWRAP_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include "jlcxx/jlcxx.hpp"

#pragma GCC diagnostic pop

JLCXX_MODULE define_julia_module(jlcxx::Module& mod);

#endif   // CXXWRAP_H
