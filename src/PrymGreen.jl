module PrymGreen

using LightXML
using Cxx
using Singular

export load_example, submatrix, check_prym_green_conjecture

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
    R, X = PolynomialRing(Fp(char), vars; degree_bound = 3)
    global X
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    I = Ideal(R, [ eval(parse(s)) for s in basis ])
    R, I, key
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
    ordstr = icxx"""rOrdStr($ring);"""
println(unsafe_string(ordstr));
    if !ismatch(r"^dp\([0-9].*\),c", unsafe_string(ordstr))
        error("monomial ordering must be (dp, c)")
    end
    icxx"""omFree($ordstr);"""
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

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
