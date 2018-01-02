oldwdir = pwd()

pkgdir = dirname(dirname(@__FILE__))
wdir = joinpath(pkgdir, "deps")
vdir = joinpath(pkgdir, "local")
sdir = joinpath(pkgdir, "src")
ldir = joinpath(sdir, "lib")

cd(ldir)

run(`./configure --prefix=$vdir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
