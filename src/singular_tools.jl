using Cxx
using Singular
using Singular.libSingular

function rOrdStr(r::libSingular.ring)
    ordstr = icxx"""rOrdStr($r);"""
    res = unsafe_string(ordstr)
    icxx"""omFree($ordstr);"""
    res
end
