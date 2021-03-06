@info("Building PrymGreen")

pkg_dir = realpath(joinpath(@__DIR__, ".."))
lib_dir = joinpath(pkg_dir, "src", "lib")
local_dir = joinpath(pkg_dir, "local")
cxxwrap_dir = Base.find_package("PrymGreen", "CxxWrap")
cxxwrap_dir = realpath(joinpath(dirname(cxxwrap_dir), "..", "deps", "usr"))
flint_dir = Base.find_package("PrymGreen", "Nemo")
flint_dir = realpath(joinpath(dirname(flint_dir), "..", "deps", "usr"))
julia_dir = realpath(joinpath(Sys.BINDIR, ".."))
singular_dir = Base.find_package("PrymGreen", "Singular")
singular_dir = realpath(joinpath(dirname(singular_dir), "..", "local"))

oldwdir = pwd()
cd(lib_dir)

run(`./configure --enable-silent-rules --prefix=$local_dir \
    --with-cxxwrap=$cxxwrap_dir \
    --with-flint=$flint_dir \
    --with-julia=$julia_dir \
    --with-singular=$singular_dir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
