# PrymGreen

[![Build Status](https://travis-ci.com/steenpass/PrymGreen.jl.svg?branch=master)](https://travis-ci.com/steenpass/PrymGreen.jl)
[![Coverage Status](https://coveralls.io/repos/github/steenpass/PrymGreen.jl/badge.svg)](https://coveralls.io/github/steenpass/PrymGreen.jl)
[![codecov](https://codecov.io/gh/steenpass/PrymGreen.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/steenpass/PrymGreen.jl)

A Julia package for computations related to the Prym-Green conjecture.

To build PrymGreen.jl, start Julia and then type:

```julia
julia> using Pkg
julia> singular_url = "https://github.com/oscar-system/Singular.jl"
julia> Pkg.add(PackageSpec(url = singular_url))
julia> prymgreen_url = "https://github.com/steenpass/PrymGreen.jl"
julia> Pkg.add(PackageSpec(url = prymgreen_url))
```
To use PrymGreen.jl, start Julia and then type:
```julia
julia> using PrymGreen
```
