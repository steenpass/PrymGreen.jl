oldwdir = pwd()

pkgdir = dirname(dirname(@__FILE__))
wdir = joinpath(pkgdir, "deps")
vdir = joinpath(pkgdir, "local")
sdir = joinpath(pkgdir, "src")
ldir = joinpath(sdir, "lib")
flint_dir = joinpath(Pkg.dir("Nemo"), "local", "lib")

cd(ldir)

run(`./configure --prefix=$vdir --with-flint=$flint_dir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
