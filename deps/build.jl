@info("Building PrymGreen")

import Nemo
import Singular

pkg_dir = realpath(joinpath(@__DIR__, ".."))
lib_dir = joinpath(pkg_dir, "src", "lib")
local_dir = joinpath(pkg_dir, "local")
flint_dir = realpath(joinpath(dirname(pathof(Nemo)), "..", "local"))
singular_dir = realpath(joinpath(dirname(pathof(Singular)), "..", "local"))

oldwdir = pwd()
cd(lib_dir)

run(`./configure --enable-silent-rules --prefix=$local_dir \
    --with-flint=$flint_dir --with-singular=$singular_dir`)
run(`make`)
run(`make install`)
run(`make clean`)

cd(oldwdir)
