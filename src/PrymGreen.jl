module PrymGreen

using LightXML
using Cxx
using Singular

export load_example, resort, submatrix, run_example,
    check_prym_green_conjecture

const pkgdir = dirname(dirname(@__FILE__))
addHeaderDir(joinpath(pkgdir, "local", "include"), kind = C_User)
cxxinclude(joinpath("submatrix.h"), isAngled = false)

global const Msize_t = typeof(icxx"""(msize_t)0;""")
global const Nvals_t = typeof(icxx"""(nvals_t)0;""")
global const Entry_t = typeof(icxx"""(entry_t)0;""")

function __init__()
    ldir = joinpath(pkgdir, "local", "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
    Libdl.dlopen(joinpath("libprymgreen"), Libdl.RTLD_GLOBAL)
end

include("singular_tools.jl")

function load_example(filename::String)
    file = open(filename)
    s = "<root>\n" * readstring(file) * "</root>"
    close(file)
    example_xml = parse_string(s)
    root_xml = root(example_xml)
    key = LightXML.content(find_element(root_xml, "Key"))
    char = parse(Int64, LightXML.content(find_element(root_xml, "char")))
    vars = LightXML.content(find_element(root_xml, "vars"))
    vars = split(vars, ['\n', ' ', '[', ',', ']'], keep=false)
    vars = Array{String,1}(vars)
    basis = LightXML.content(find_element(root_xml, "basis"))
    basis = split(basis, ['\n', ' ', '[', ',', ']'], keep=false)
    basis = Array{String,1}(basis)
    free(example_xml)
    global X
    R, X = PolynomialRing(Fp(char), vars)
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    I = Ideal(R, [ eval(parse(s)) for s in basis ])
    R, I, key
end

function resort(I::sideal)
    R = Singular.base_ring(I)
    J = Ideal(R, [R() for i in 1:ngens(I)]...)
    i = 1
    for j = 1:ngens(I)
        if degree(I[j]) != degree(I[i])
            for k = (j-1):-1:i
                J[i-1+j-k] = I[k]
            end
            i = j
        end
    end
    for k = ngens(I):-1:i
        J[i+ngens(I)-k] = I[k]
    end
    J
end

#=
This function can be simplified as soon as maps between polynomial rings are
available in Singular.jl.
=#
function set_degree_bound(R::PolyRing, I::sideal, d::Int)
    for i in 1:ngens(I)
        if degree(I[i]) > d
            error("degree bound too low")
        end
    end
    I.isGB = true
    r = fres(I, 0, "frame")
    B = betti(r)
    if size(B, 1) > d
        error("degree bound possibly too low to compute free resolution")
    end
    vars = [ string(Singular.gens(R)[i]) for i in 1:ngens(R) ]
    global X
    S, X = PolynomialRing(base_ring(R), vars; degree_bound = d)
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    Ideal(S, [ eval(parse(string(I[i]))) for i in 1:ngens(I) ])
end

function submatrix(r::Singular.sresolution, g::Int, char::Entry_t)
    index = div(g, 2)-2
    B = betti(r)
    size = Msize_t(B[3, index])
    if size != Msize_t(B[2, index+1])
        error("matrix not square")
    end
println("size = ", size)
    limit = Msize_t(B[2, index])
    r_ptr = r.ptr
    R = Singular.base_ring(r)
    ring = R.ptr
    ordstr = rOrdStr(ring)
println(ordstr);
    if !ismatch(r"^dp\([0-9].*\),c", ordstr)
        error("monomial ordering must be (dp, c)")
    end
    values_ptr = icxx"""(entry_t **)malloc(sizeof(entry_t *));"""
    n_values = icxx"""
            check_matrix($values_ptr, $r_ptr, $g, $size, $limit, $char, $ring);
            """
println("n_values = ", n_values)
    if n_values == 0
        error("number of values in prym green matrix must be positive")
    end
    A = unsafe_wrap(Array, unsafe_load(values_ptr), (n_values, ), true)
    icxx"""free($values_ptr);"""
    A
end

function run_example(filename::String)
    R, I, key = load_example(filename)
println(key)
    g = parse(Int, match(r"(?<=g)(.*)(?=_)", key).match)
println("g = ", g)
    char = parse(PrymGreen.Entry_t, match(r"(?<=@)(.*)(?=g)", key).match)
println("char = ", char)
    I = std(I; complete_reduction = true)
    I = resort(I)
    I = set_degree_bound(R, I, 3)
    I.isGB = true
    gc()
    @time r = fres(I, div(g, 2)-2, "single module")
println(r)
    @time A = submatrix(r, g, char)
println(map(x -> Int(x), A[1:10]))
end

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
