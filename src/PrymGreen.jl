module PrymGreen

export mytest

function __init__()
    pkgdir = dirname(dirname(@__FILE__))
    vdir = joinpath(pkgdir, "local")
    ldir = joinpath(vdir, "lib")
    push!(Libdl.DL_LOAD_PATH, ldir)
end

function mytest(x, y)
    return x + y
end

end # module
