using Cxx
using Nemo
using Singular
using Singular.libSingular

#=
modified copy from Singular.jl (commit b277b7c) to make use_cache and
use_tensor_trick available
=#
function id_fres(I::libSingular.ideal, n::Cint, method::String,
        R::libSingular.ring, use_cache::Bool, use_tensor_trick::Bool)
    s = icxx"""const ring origin = currRing;
            rChangeCurrRing($R);
            syStrategy s = syFrank($I, $n, $method, $use_cache,
                    $use_tensor_trick);
            rChangeCurrRing(origin);
            s;
        """
    r = icxx"""$s->fullres;"""
    length = icxx"""$s->length;"""
    r, Int(length)
end

#=
modified copy from Singular.jl (commit b277b7c) to make use_cache and
use_tensor_trick available
=#
function fres{T <: Nemo.RingElem}(id::sideal{T}, max_length::Int,
        method::String = "complete";
        use_cache::Bool = true, use_tensor_trick::Bool = false)
    id.isGB == false && error("ideal is not a standard basis")
    max_length < 0 && error("length for fres must not be negative")
    R = base_ring(id)
    if max_length == 0
        max_length = ngens(R)+1
        # TODO: consider qrings
    end
    if (method != "complete"
            && method != "frame"
            && method != "extended frame"
            && method != "single module")
        error("wrong optional argument for fres")
    end
    r, length = id_fres(id.ptr, Cint(max_length), method, R.ptr, use_cache,
            use_tensor_trick)
    return Singular.sresolution{T}(R, length, r)
end

function rOrdStr(r::libSingular.ring)
    ordstr = icxx"""rOrdStr($r);"""
    res = unsafe_string(ordstr)
    icxx"""omFree($ordstr);"""
    res
end
