module PrymGreen

import AbstractAlgebra
import CxxWrap
import DelimitedFiles
import Hecke
import Libdl
import LightXML
import Nemo
import Random
import Singular

export init_rng, load_example, test_example, run_example,
        test_modular_arithmetic

const pkgdir = realpath(joinpath(@__DIR__, ".."))
const libdir = realpath(joinpath(pkgdir, "local", "lib"))

CxxWrap.@wrapmodule(realpath(joinpath(libdir, "libprymgreen.so")))

function __init__()
    ldir = joinpath(pkgdir, "local", "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
    Libdl.dlopen("libprymgreen.so", Libdl.RTLD_GLOBAL)
    CxxWrap.@initcxx
end

global const Msize_t = typeof(return_msize_t())
global const Nvals_t = typeof(return_nvals_t())
global const Entry_t = typeof(return_entry_t())
global const Arith_t = typeof(return_arith_t())

include("singular_tools.jl")
include("generate_random_pcnc.jl")

macro time_info(ex)
    return :($(esc(:(print_info))) ? $(esc(:(@time $ex))) : $(esc(ex)))
end

function glibc(S::Symbol)
    libname = "libc.so.6"
    @static if Sys.isapple()
        libname = "libSystem"
    end
    return Libdl.dlsym(Libdl.dlopen(libname), S)
end

function malloc(::Type{T}, number::Unsigned = UInt(1)) where T
    return ccall(glibc(:malloc), Ptr{T}, (Csize_t, ), number*sizeof(T))
end

function free(object::T) where T
    ccall(glibc(:free), Nothing, (T, ), object)
end

function read_xml_file(filename::String)
    file = open(filename)
    s = "<root>\n" * read(file, String) * "</root>"
    close(file)
    example_xml = LightXML.parse_string(s)
    root_xml = LightXML.root(example_xml)
    key = LightXML.content(LightXML.find_element(root_xml, "Key"))
    char = LightXML.content(LightXML.find_element(root_xml, "char"))
    char = parse(Entry_t, char)
    vars = LightXML.content(LightXML.find_element(root_xml, "vars"))
    vars = split(vars, ['\n', ' ', '[', ',', ']'], keepempty = false)
    vars = Array{String, 1}(vars)
    basis = LightXML.content(LightXML.find_element(root_xml, "basis"))
    basis = split(basis, ['\n', ' ', '[', ',', ']'], keepempty = false)
    basis = Array{String, 1}(basis)
    LightXML.free(example_xml)
    return (key, char, vars, basis)
end

function load_example(filename::String, print_info::Bool = false)
    key, char, vars, basis = read_xml_file(filename)
    print_info && println(key)
    global X
    R, X = Singular.PolynomialRing(Singular.Fp(Int(char)), vars;
            ordering = :degrevlex, ordering2 = :comp1max)
    for (i, s) in enumerate(vars)
        eval(Meta.parse("$s = X[$i]"))
    end
    I = Singular.Ideal(R, [ eval(Meta.parse(s)) for s in basis ])
    I = Singular.std(I; complete_reduction = true)
    I = resort(I)
    g = parse(Int, match(r"(?<=g)(.*)(?=_)", key).match)
    char = parse(PrymGreen.Entry_t, match(r"(?<=@)(.*)(?=g)", key).match)
    return (R, I, char, g)
end

function rev_dp(a::Singular.spoly{T}, b::Singular.spoly{T}) where
        T <: AbstractAlgebra.RingElem
    deg_a = Singular.total_degree(a)
    deg_b = Singular.total_degree(b)
    if (deg_a == deg_b)
        R = Singular.parent(a)
        return (p_LmCmp(a.ptr, b.ptr, R.ptr) == 1)
    end
    return (deg_a < deg_b)
end

function resort(I::Singular.sideal{T}) where T <: AbstractAlgebra.RingElem
    R = Singular.base_ring(I)
    gens = [ I[i] for i = 1:Singular.ngens(I) ]
    J = Singular.Ideal(R, sort!(gens; lt=rev_dp))
    J.isGB = I.isGB
    return J
end

#=
This function can be simplified as soon as maps between polynomial rings are
available in Singular.jl.
=#
function set_degree_bound(R::Singular.PolyRing, I::Singular.sideal, d::Int)
    for i in 1:Singular.ngens(I)
        if Singular.total_degree(I[i]) > d
            error("degree bound too low")
        end
    end
    r = PrymGreen.fres(I, 0, "frame")
    B = Singular.betti(r)
    if size(B, 1) > d
        error("degree bound possibly too low to compute free resolution")
    end
    vars = [ string(Singular.gens(R)[i]) for i in 1:Singular.nvars(R) ]
    global X
    S, X = Singular.PolynomialRing(Singular.base_ring(R), vars;
            degree_bound = d, ordering = :degrevlex, ordering2 = :comp1max)
    for (i, s) in enumerate(vars)
        eval(Meta.parse("$s = X[$i]"))
    end
    J = Singular.Ideal(S,
            [ eval(Meta.parse(string(I[i]))) for i in 1:Singular.ngens(I) ])
    J.isGB = I.isGB
    return (S, J)
end

function check_ordering(R::Singular.PolyRing)
    ordstr = rOrdStr(R.ptr)
    if !occursin(r"^dp\([0-9].*\),c", ordstr)
        error("monomial ordering must be (dp, c)")
    end
end

function betti_table_entries(res::Singular.sresolution, g::Int)
    index = div(g, 2)-2
    B = Singular.betti(res)
    prym_green_size = Msize_t(B[3, index])
    if prym_green_size != Msize_t(B[2, index+1])
        error("matrix not square")
    end
    limit = Msize_t(B[2, index])
    return (prym_green_size, limit)
end

function submatrix(res::Singular.sresolution, R::Singular.PolyRing, g::Int,
        char::Entry_t)
    check_ordering(R)
    prym_green_size, limit = betti_table_entries(res, g)
    values_ptr = malloc(Ptr{Entry_t})
    n_values = check_matrix(Ptr{Nothing}(values_ptr), res.ptr, g,
            prym_green_size, limit, char, R.ptr)
    if n_values == 0
        error("number of values in Prym-Green matrix must be positive")
    end
    A = unsafe_wrap(Array{Entry_t, 1}, unsafe_load(values_ptr), (n_values, );
            own = true)
    free(values_ptr)
    return (A, prym_green_size)
end

function init_rng()
    seed = rand(Random.RandomDevice(), UInt32, 4)
    println("rng seed = ", seed)
    return Random.MersenneTwister(seed)
end

function multiply_matrix(A::Array{Arith_t, 1}, v::Array{Arith_t, 1}, g::Int,
        char::Entry_t)
    char = Arith_t(char)
    Axv_ptr = malloc(Ptr{Arith_t})
    length_Axv = ccall((:multiply_matrix, "libprymgreen.so"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Ptr{Arith_t}, Int, Arith_t),
            Axv_ptr, A, v, g, char)
    Axv = unsafe_wrap(Array{Arith_t, 1}, unsafe_load(Axv_ptr), (length_Axv, );
            own = true)
    free(Axv_ptr)
    return Axv
end

function dense_pg_matrix(res::Singular.sresolution, R::Singular.PolyRing,
        g::Int, prym_green_size::Msize_t, limit::Msize_t)
    A_dense_ptr = malloc(Ptr{Entry_t}, prym_green_size)
    size_A = dense_matrix(Ptr{Nothing}(A_dense_ptr), res.ptr, g,
            prym_green_size, limit, R.ptr)
    A_dense = unsafe_wrap(Array{Entry_t, 2}, unsafe_load(A_dense_ptr),
            (size_A, size_A); own = true)
    free(A_dense_ptr)
    return A_dense
end

function gauss(A_dense::Array{Entry_t, 2}, char::Entry_t)
    @time rank = ccall((:gauss, "libprymgreen.so"), Msize_t,
            (Ptr{Entry_t}, Msize_t, Msize_t, Entry_t),
            A_dense, size(A_dense, 1), size(A_dense, 2), char)
    println("rank:      ", rank)
end

function write_dense_matrix(A_dense::Array{Entry_t, 2}, g::Int, char::Entry_t)
    file = "pcnc_g"*string(g)*"_submatrix_"*string(size(A_dense, 1))
    if isfile(file)
        rm(file)
    end
    f = open(file, "a")
    write(f, string(char)*"\n")
    write(f, string(size(A_dense, 1))*" "*string(size(A_dense, 2))*"\n")
    A_print = Array{Int, 2}(A_dense)
    DelimitedFiles.writedlm(f, A_print, ' ')
    close(f)
end

function check_multiplication(A::Array{Entry_t, 1}, res::Singular.sresolution,
        R::Singular.PolyRing, g::Int, char::Entry_t, rng::Random.AbstractRNG,
        print_info::Bool = false)
    g > 18 && return
    prym_green_size, limit = betti_table_entries(res, g)
    A = Array{Arith_t, 1}(A)
    v = Array{Arith_t, 1}(rand(rng, 0:(char-1), Int(prym_green_size)))
    Axv = multiply_matrix(A, v, g, char)
    A_dense = dense_pg_matrix(res, R, g, prym_green_size, limit)
    g <= 12 && gauss(A_dense, char)
    g <= 12 && write_dense_matrix(A_dense, g, char)
    A_dense = Array{UInt128, 2}(A_dense)
    print_info && print("mlt. test: ")
    if Axv == A_dense*v .% char
        print_info && printstyled("passed\n"; color = :green)
        return true
    else
        print_info && printstyled("failed\n"; bold = true, color = :red)
        return false
    end
end

function print_matrix_info(A::Array{Entry_t, 1}, prym_green_size::Msize_t)
    println("p_g_size = ", prym_green_size)
    println("n_values = ", size(A, 1))
    println("max bytes: ", Singular_MaxBytesSystem())
end

function recurrence_sequence(A::Array{Entry_t, 1}, prym_green_size::Msize_t,
        g::Int, char::Entry_t, rng::Random.AbstractRNG)
    A = Array{Arith_t, 1}(A)
    char = Arith_t(char)
    v = Array{Arith_t, 1}(rand(rng, 0:(char-1), Int(prym_green_size)))
    index = Msize_t(rand(rng, 0:(prym_green_size-1)))
    seq_ptr = malloc(Ptr{Arith_t})
    length_seq = ccall((:recurrence_sequence, "libprymgreen.so"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Nvals_t, Ptr{Arith_t}, Msize_t,
                Msize_t, Int, Arith_t),
            seq_ptr, A, size(A, 1), v, prym_green_size, index, g, char)
    seq = unsafe_wrap(Array{Arith_t, 1}, unsafe_load(seq_ptr), (length_seq, );
            own = true)
    free(seq_ptr)
    return seq
end

function berlekamp_massey(S::Array{Arith_t, 1}, char::Entry_t)
    char = Arith_t(char)
    lfsr_ptr = malloc(Ptr{Arith_t})
    length_lfsr = ccall((:berlekamp_massey, "libprymgreen.so"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Msize_t, Arith_t),
            lfsr_ptr, S, size(S, 1), char)
    lfsr = unsafe_wrap(Array{Arith_t, 1}, unsafe_load(lfsr_ptr),
            (length_lfsr, ); own = true)
    free(lfsr_ptr)
    return lfsr
end

function check_berlekamp_massey(C::Array{Arith_t, 1}, S::Array{Arith_t, 1},
        char::Entry_t, print_info::Bool = false)
    R = Nemo.NmodRing(UInt64(char))
    Sp = [ R(i) for i in S ]
    f = Hecke.berlekamp_massey(Sp)
    C_check = [ Arith_t(Nemo.coeff(f, i).data) for i in Nemo.degree(f):-1:0 ]
    print_info && print("B-M. test: ")
    if C == C_check
        print_info && printstyled("passed\n"; color = :green)
        return true
    else
        print_info && printstyled("failed\n"; bold = true, color = :red)
        return false
    end
end

function kernel(C::Array{Arith_t, 1}, A::Array{Entry_t, 1},
        prym_green_size::Msize_t, g::Int, char::Entry_t,
        rng::Random.AbstractRNG)
    A = Array{Arith_t, 1}(A)
    char = Arith_t(char)
    v = Array{Arith_t, 1}(rand(rng, 0:(char-1), Int(prym_green_size)))
    ker_ptr = malloc(Ptr{Arith_t})
    length_ker = ccall((:kernel, "libprymgreen.so"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Ptr{Arith_t}, Nvals_t,
                Ptr{Arith_t}, Msize_t, Int, Arith_t),
            ker_ptr, C, A, size(A, 1), v, prym_green_size, g, char)
    ker = unsafe_wrap(Array{Arith_t, 1}, unsafe_load(ker_ptr), (length_ker, );
            own = true)
    free(ker_ptr)
    return ker
end

function test_example(R::Singular.PolyRing, I::Singular.sideal, char::Entry_t,
        g::Int, rng::Random.AbstractRNG; print_info::Bool = false)
    success = true
    R, I = set_degree_bound(R, I, 3)
    GC.gc()
    @time_info res = PrymGreen.fres(I, div(g, 2)-2, "single module";
            use_cache = false, use_tensor_trick = true)
    @time_info A, prym_green_size = submatrix(res, R, g, char)
    print_info && print_matrix_info(A, prym_green_size)
    print_info && println("A[1:4]   = ", map(x -> Int(x), A[1:4]))
    success &= check_multiplication(A, res, R, g, char, rng, print_info)
    res = nothing
    GC.gc()
    @time_info S = recurrence_sequence(A, prym_green_size, g, char, rng)
    print_info && println("S[1:4]   = ", map(x -> Int(x), S[1:4]))
    @time_info C = PrymGreen.berlekamp_massey(S, char)
    print_info && println("C[1:4]   = ", map(x -> Int(x), C[1:4]))
    success &= check_berlekamp_massey(C, S, char, print_info)
    return success
end

function run_example(R::Singular.PolyRing, I::Singular.sideal, char::Entry_t,
        g::Int, rng::Random.AbstractRNG; print_info::Bool = false)
    R, I = set_degree_bound(R, I, 3)
    GC.gc()
    @time_info res = PrymGreen.fres(I, div(g, 2)-2, "single module";
            use_cache = false, use_tensor_trick = true)
    @time_info A, prym_green_size = submatrix(res, R, g, char)
    res = nothing
    GC.gc()
    print_info && print_matrix_info(A, prym_green_size)
    @time_info S = recurrence_sequence(A, prym_green_size, g, char, rng)
    @time_info C = PrymGreen.berlekamp_massey(S, char)
    if C[end] == 0
        return false
    end
    if size(C, 1)-1 == prym_green_size
        return true
    end
    return nothing
end

function test_modular_arithmetic()
    @time res = ccall((:mult_preinv_test, "libprymgreen.so"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 100000000)
    return res
end

end # module
