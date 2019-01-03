module PrymGreen

import AbstractAlgebra
import CxxWrap
import Libdl
import LightXML

import Nemo
import Hecke
import Singular

export run_example, check_prym_green_conjecture

const pkgdir = realpath(joinpath(@__DIR__, ".."))
const libdir = realpath(joinpath(pkgdir, "local", "lib"))

CxxWrap.@wrapmodule(realpath(joinpath(libdir, "libprymgreen." * Libdl.dlext)))

function __init__()
    CxxWrap.@initcxx
end

global const Msize_t = typeof(return_msize_t())
global const Nvals_t = typeof(return_nvals_t())
global const Entry_t = typeof(return_entry_t())
global const Arith_t = typeof(return_arith_t())

include("singular_tools.jl")
include("generate_random_pcnc.jl")

function read_xml_file(filename::String)
    file = open(filename)
    s = "<root>\n" * readstring(file) * "</root>"
    close(file)
    example_xml = LightXML.parse_string(s)
    root_xml = LightXML.root(example_xml)
    key = LightXML.content(LightXML.find_element(root_xml, "Key"))
    char = LightXML.content(LightXML.find_element(root_xml, "char"))
    char = parse(Entry_t, char)
    vars = LightXML.content(LightXML.find_element(root_xml, "vars"))
    vars = split(vars, ['\n', ' ', '[', ',', ']'], keep = false)
    vars = Array{String, 1}(vars)
    basis = LightXML.content(LightXML.find_element(root_xml, "basis"))
    basis = split(basis, ['\n', ' ', '[', ',', ']'], keep = false)
    basis = Array{String, 1}(basis)
    LightXML.free(example_xml)
    key, char, vars, basis
end

function load_example(filename::String)
    key, char, vars, basis = read_xml_file(filename)
    global X
    R, X = Singular.PolynomialRing(Singular.Fp(Int(char)), vars;
            ordering = :degrevlex, ordering2 = :comp1max)
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    I = Singular.Ideal(R, [ eval(parse(s)) for s in basis ])
    R, I, key
end

function resort(I::Singular.sideal, R::Singular.PolyRing)
    J = Singular.Ideal(R, [ R() for i in 1:Singular.ngens(I) ]...)
    i = 1
    for j = 1:Singular.ngens(I)
        if Singular.degree(I[j]) != Singular.degree(I[i])
            for k = (j-1):-1:i
                J[i+j-1-k] = I[k]
            end
            i = j
        end
    end
    for k = Singular.ngens(I):-1:i
        J[i+Singular.ngens(I)-k] = I[k]
    end
    J
end

#=
This function can be simplified as soon as maps between polynomial rings are
available in Singular.jl.
=#
function set_degree_bound(R::Singular.PolyRing, I::Singular.sideal, d::Int)
    for i in 1:Singular.ngens(I)
        if Singular.degree(I[i]) > d
            error("degree bound too low")
        end
    end
    I.isGB = true
    r = PrymGreen.fres(I, 0, "frame")
    B = Singular.betti(r)
    if size(B, 1) > d
        error("degree bound possibly too low to compute free resolution")
    end
    vars = [ string(Singular.gens(R)[i]) for i in 1:Singular.ngens(R) ]
    global X
    S, X = Singular.PolynomialRing(Singular.base_ring(R), vars;
            degree_bound = d, ordering = :degrevlex, ordering2 = :comp1max)
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    J = Singular.Ideal(S,
            [ eval(parse(string(I[i]))) for i in 1:Singular.ngens(I) ])
    S, J
end

function check_ordering(R::Singular.PolyRing)
    ordstr = rOrdStr(R.ptr)
    if !ismatch(r"^dp\([0-9].*\),c", ordstr)
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
    prym_green_size, limit
end

function submatrix(res::Singular.sresolution, R::Singular.PolyRing, g::Int,
        char::Entry_t)
    check_ordering(R)
    prym_green_size, limit = betti_table_entries(res, g)
    values_ptr = ccall((:malloc, "libc"), Ptr{Ptr{Entry_t}}, (Csize_t, ),
            sizeof(Ptr{Entry_t}))
    n_values = check_matrix(values_ptr, res.ptr, g, prym_green_size, limit,
        char, R.ptr)
    if n_values == 0
        error("number of values in Prym-Green matrix must be positive")
    end
    A = unsafe_wrap(Array, unsafe_load(values_ptr), (n_values, ), true)
    ccall((:free, "libc"), Nothing, (Ptr{Ptr{Entry_t}}, ), values_ptr)
    A, prym_green_size
end

function init_rng()
    seed = rand(RandomDevice(), UInt32, 4)
    println("rng seed = ", seed)
    MersenneTwister(seed)
end

function multiply_matrix(A::Array{Arith_t, 1}, v::Array{Arith_t, 1}, g::Int,
        char::Entry_t)
    char = Arith_t(char)
    Axv_ptr = ccall((:malloc, "libc"), Ptr{Ptr{Arith_t}}, (Csize_t, ),
            sizeof(Ptr{Arith_t}))
    length_Axv = ccall((:multiply_matrix, "libprymgreen"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Ptr{Arith_t}, Int, Arith_t),
            Axv_ptr, A, v, g, char)
    Axv = unsafe_wrap(Array, unsafe_load(Axv_ptr), (length_Axv, ), true)
    ccall((:free, "libc"), Nothing, (Ptr{Ptr{Arith_t}}, ), Axv_ptr)
    Axv
end

function dense_pg_matrix(res::Singular.sresolution, R::Singular.PolyRing,
        g::Int, prym_green_size::Msize_t, limit::Msize_t)
    A_dense_ptr = ccall((:malloc, "libc"), Ptr{Ptr{Entry_t}}, (Csize_t, ),
            prym_green_size*sizeof(Ptr{Entry_t}))
    size_A = dense_matrix(A_dense_ptr, res.ptr, g, prym_green_size, limit,
            R.ptr);
    A_dense = unsafe_wrap(Array, unsafe_load(A_dense_ptr), (size_A, size_A),
            true)
    ccall((:free, "libc"), Nothing, (Ptr{Ptr{Entry_t}}, ), A_dense_ptr)
    A_dense
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
    writedlm(f, A_print, ' ')
    close(f)
end

function gauss(A_dense::Array{Entry_t, 2}, char::Entry_t)
    @time rank = ccall((:gauss, "libprymgreen"), Msize_t,
            (Ptr{Entry_t}, Msize_t, Msize_t, Entry_t),
            A_dense, size(A_dense, 1), size(A_dense, 2), char)
    println("rank:      ", rank)
end

function check_multiplication(A::Array{Entry_t, 1}, res::Singular.sresolution,
        R::Singular.PolyRing, g::Int, char::Entry_t, rng::AbstractRNG)
    g > 18 && return
    prym_green_size, limit = betti_table_entries(res, g)
    A = Array{Arith_t, 1}(A)
    v = Array{Arith_t, 1}(rand(rng, 0:(char-1), Int(prym_green_size)))
    Axv = multiply_matrix(A, v, g, char)
    A_dense = dense_pg_matrix(res, R, g, prym_green_size, limit)
    A_dense = Array{UInt128, 2}(A_dense)
    # gauss(A_dense, char)
    # write_dense_matrix(A_dense, g, char)
    print("mlt. test: ")
    if Axv == A_dense*v .% char
        print_with_color(:green, "passed\n")
    else
        print_with_color(:red, "failed\n"; bold = true)
    end
end

function print_matrix_info(A::Array{Entry_t, 1}, prym_green_size::Msize_t)
    println("p_g_size = ", prym_green_size)
    println("n_values = ", size(A, 1))
    println("max bytes: ", Singular_MaxBytesSystem())
    println("A[1:4]   = ", map(x -> Int(x), A[1:4]))
end

function recurrence_sequence(A::Array{Entry_t, 1}, prym_green_size::Msize_t,
        g::Int, char::Entry_t, rng::AbstractRNG)
    A = Array{Arith_t, 1}(A)
    char = Arith_t(char)
    v = Array{Arith_t, 1}(rand(rng, 0:(char-1), Int(prym_green_size)))
    index = Msize_t(rand(rng, 0:(prym_green_size-1)))
    seq_ptr = ccall((:malloc, "libc"), Ptr{Ptr{Arith_t}}, (Csize_t, ),
            sizeof(Ptr{Arith_t}))
    length_seq = ccall((:recurrence_sequence, "libprymgreen"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Nvals_t, Ptr{Arith_t}, Msize_t,
                Msize_t, Int, Arith_t),
            seq_ptr, A, size(A, 1), v, prym_green_size, index, g, char)
    seq = unsafe_wrap(Array, unsafe_load(seq_ptr), (length_seq, ), true)
    ccall((:free, "libc"), Void, (Ptr{Ptr{Arith_t}}, ), seq_ptr)
    seq
end

function berlekamp_massey(S::Array{Arith_t, 1}, char::Entry_t)
    char = Arith_t(char)
    lfsr_ptr = ccall((:malloc, "libc"), Ptr{Ptr{Arith_t}}, (Csize_t, ),
            sizeof(Ptr{Arith_t}))
    length_lfsr = ccall((:berlekamp_massey, "libprymgreen"), Msize_t,
            (Ptr{Ptr{Arith_t}}, Ptr{Arith_t}, Msize_t, Arith_t),
            lfsr_ptr, S, size(S, 1), char)
    lfsr = unsafe_wrap(Array, unsafe_load(lfsr_ptr), (length_lfsr, ), true)
    ccall((:free, "libc"), Void, (Ptr{Ptr{Arith_t}}, ), lfsr_ptr)
    lfsr
end

function check_berlekamp_massey(C::Array{Arith_t, 1}, S::Array{Arith_t, 1},
        char::Entry_t)
    R = Nemo.NmodRing(UInt64(char))
    Sp = [ R(i) for i in S ]
    f = Hecke.berlekamp_massey(Sp)
    C_check = [ Arith_t(Nemo.coeff(f, i).data) for i in Nemo.degree(f):-1:0 ]
    print("B-M. test: ")
    if C == C_check
        print_with_color(:green, "passed\n")
    else
        print_with_color(:red, "failed\n"; bold = true)
    end
end

function run_example(filename::String; print_info::Bool = false)
    R, I, key = load_example(filename)
    print_info && println(key)
    g = parse(Int, match(r"(?<=g)(.*)(?=_)", key).match)
    char = parse(PrymGreen.Entry_t, match(r"(?<=@)(.*)(?=g)", key).match)
    I = Singular.std(I; complete_reduction = true)
    I = resort(I, R)
    R, I = set_degree_bound(R, I, 3)
    I.isGB = true
    gc()
    @time res = PrymGreen.fres(I, div(g, 2)-2, "single module";
            use_cache = false, use_tensor_trick = true)
    @time A, prym_green_size = submatrix(res, R, g, char)
    print_info && print_matrix_info(A, prym_green_size)
    rng = init_rng()
    check_multiplication(A, res, R, g, char, rng)
    res = nothing
    gc()
    @time S = recurrence_sequence(A, prym_green_size, g, char, rng)
    print_info && println("S[1:4]   = ", map(x -> Int(x), S[1:4]))
    @time C = PrymGreen.berlekamp_massey(S, char)
    print_info && println("C[1:4]   = ", map(x -> Int(x), C[1:4]))
    check_berlekamp_massey(C, S, char)
    nothing
end

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
