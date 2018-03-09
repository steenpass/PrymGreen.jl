info("Building PrymGreen")

pkg_dir = Pkg.dir("PrymGreen")
lib_dir = joinpath(pkg_dir, "src", "lib")
local_dir = joinpath(pkg_dir, "local")
flint_dir = joinpath(Pkg.dir("Nemo"), "local")
singular_dir = joinpath(Pkg.dir("Singular"), "local")

oldwdir = pwd()
cd(lib_dir)

run(`./configure --enable-silent-rules --prefix=$local_dir \
    --with-flint=$flint_dir --with-singular=$singular_dir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
