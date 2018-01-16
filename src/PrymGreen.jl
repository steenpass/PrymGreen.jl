module PrymGreen

using LightXML
using Singular

export load_example, check_prym_green_conjecture

function __init__()
    pkgdir = dirname(dirname(@__FILE__))
    ldir = joinpath(pkgdir, "local", "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
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

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
