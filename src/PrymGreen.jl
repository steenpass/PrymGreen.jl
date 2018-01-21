module PrymGreen

using LightXML
using Cxx
using Singular

export load_example, submatrix, prym_green_matrix, check_prym_green_conjecture

function __init__()
    pkgdir = dirname(dirname(@__FILE__))
    ldir = joinpath(pkgdir, "local", "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
    Libdl.dlopen(joinpath("libprymgreen"), Libdl.RTLD_GLOBAL)
    addHeaderDir(joinpath(pkgdir, "local", "include"), kind = C_User)
    cxxinclude(joinpath("submatrix.h"), isAngled = false)
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
    R, X = PolynomialRing(Fp(char), vars)
    global X
    for (i, s) in enumerate(vars)
        eval(parse("$s = X[$i]"))
    end
    I = Ideal(R, [ eval(parse(s)) for s in basis ])
    R, I, key
end

function submatrix(r::Singular.sresolution, g::Int, char::Int)
    index = div(g, 2)-2
    B = betti(r)
    size = Int64(B[3, index])
    if size != B[2, index+1]
        error("matrix not square")
    end
    limit = B[2, index]
    r_ptr = r.ptr
    R = Singular.base_ring(r)
    ring = R.ptr
    ordstr = icxx"""rOrdStr($ring);"""
    if !ismatch(r"^dp\([0-9].*\),c", unsafe_string(ordstr))
        error("monomial ordering must be (dp, c)")
    end
    icxx"""omFree($ordstr);"""
    values_ptr = icxx"""(int **)malloc(sizeof(int *));"""
    n_values = icxx"""
            check_matrix($values_ptr, $r_ptr, $g, $size, $limit, $char, $ring);
            """
    if n_values <= 0
        error("number of values in prym green matrix must be positive")
    end
    A = unsafe_wrap(Array, unsafe_load(values_ptr), (n_values, ), true)
    icxx"""free($values_ptr);"""
    A
end

function prym_green_matrix(r::Singular.sresolution, g::Int)
    index = div(g, 2)-2
    B = betti(r)
    size = B[3, index]
    if size != B[2, index+1]
        error("matrix not square")
    end
    limit = B[2, index]
    m = r[index]
    PR = Singular.base_ring(m)
    CR = Singular.base_ring(PR)
    A = spzeros(size, size)
    for i = 1:size
        ptr = m[i].ptr
        while ptr != C_NULL
            j = Singular.libSingular.p_GetComp(ptr, PR.ptr)-limit
            if j > 0
                n = Singular.libSingular.pGetCoeff(ptr)
                a = CR(Singular.libSingular.n_Copy(n, CR.ptr))
                A[i, j] = Int(a)
            end
            ptr = Singular.libSingular.pNext(ptr)
        end
    end
    A
end

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
