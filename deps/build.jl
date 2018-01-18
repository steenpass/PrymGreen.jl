oldwdir = pwd()

pkgdir = dirname(dirname(@__FILE__))
wdir = joinpath(pkgdir, "deps")
vdir = joinpath(pkgdir, "local")
sdir = joinpath(pkgdir, "src")
ldir = joinpath(sdir, "lib")
flint_dir = joinpath(Pkg.dir("Nemo"), "local")
singular_dir = joinpath(Pkg.dir("Singular"), "local")

cd(ldir)

run(`./configure --prefix=$vdir --with-flint=$flint_dir \
                 --with-singular=$singular_dir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
