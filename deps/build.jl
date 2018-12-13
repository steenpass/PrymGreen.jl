@info("Building PrymGreen")

import CxxWrap

@info "Loading Nemo and Singular."
@info "STDOUT and STDERR redirected to /dev/null."
s_out = stdout
s_err = stderr
redirect_stdout(open("/dev/null", "w"))
redirect_stderr(open("/dev/null", "w"))
import Nemo
import Singular
redirect_stdout(s_out)
redirect_stderr(s_err)

pkg_dir = realpath(joinpath(@__DIR__, ".."))
lib_dir = joinpath(pkg_dir, "src", "lib")
local_dir = joinpath(pkg_dir, "local")
cxxwrap_dir = realpath(joinpath(dirname(pathof(CxxWrap)), "..", "deps", "usr"))
flint_dir = realpath(joinpath(dirname(pathof(Nemo)), "..", "local"))
julia_dir = realpath(joinpath(Sys.BINDIR, ".."))
singular_dir = realpath(joinpath(dirname(pathof(Singular)), "..", "local"))

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
