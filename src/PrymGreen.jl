module PrymGreen

export check_prym_green_conjecture

function __init__()
    pkgdir = dirname(dirname(@__FILE__))
    ldir = joinpath(pkgdir, "local", "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
end

function check_prym_green_conjecture()
    @time res = ccall((:mult_preinv_test, "libprymgreen"), UInt64,
        (UInt64, UInt64, UInt64, Int64),
        3, 11, 27, 10000000000)
    return res
end

end # module
