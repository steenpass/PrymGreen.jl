using Cxx
import Nemo
import Singular

#=
modified copy from Singular.jl (commit b277b7c) to make use_cache and
use_tensor_trick available
=#
function id_fres(I::Singular.libSingular.ideal, n::Cint, method::String,
        R::Singular.libSingular.ring, use_cache::Bool, use_tensor_trick::Bool)
    s = icxx"""
            const ring origin = currRing;
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
function fres{T <: Nemo.RingElem}(id::Singular.sideal{T}, max_length::Int,
        method::String = "complete";
        use_cache::Bool = true, use_tensor_trick::Bool = false)
    id.isGB == false && error("ideal is not a standard basis")
    max_length < 0 && error("length for fres must not be negative")
    R = Singular.base_ring(id)
    if max_length == 0
        max_length = Singular.ngens(R)+1
        # TODO: consider qrings
    end
    if (method != "complete"
            && method != "frame"
            && method != "extended frame"
            && method != "single module")
        error("wrong optional argument for fres")
    end
    r, length = PrymGreen.id_fres(id.ptr, Cint(max_length), method, R.ptr,
            use_cache, use_tensor_trick)
    return Singular.sresolution{T}(R, length, r)
end

function rOrdStr(r::Singular.libSingular.ring)
    ordstr = icxx"""rOrdStr($r);"""
    res = unsafe_string(ordstr)
    icxx"""omFree($ordstr);"""
    res
end

function Singular_MaxBytesSystem()
    icxx"""
            omUpdateInfo();
            om_Info.MaxBytesSystem;
        """
end

#=
Apply the map x -> z to p where x is assumed to be a variable.
=#
function substitute_variable(p::Singular.spoly{T}, x::Singular.spoly{T},
        z::Singular.spoly{T}) where T
    R = parent(p)
    index_var = findfirst(Singular.gens(R), x)
    @assert index_var != 0
    p_ptr = p.ptr
    z_ptr = z.ptr
    R_ptr = R.ptr
    res_ptr = icxx"""
        p_SubstPoly($p_ptr, $index_var, $z_ptr, $R_ptr, $R_ptr, ndCopyMap);
    """
    return R(res_ptr)
end

#=
Successively apply the maps X[i] -> Z[i] to p where the X[i] are assumed to be
variables.
=#
function poly_substitute(p::Singular.spoly{T}, X::Array{Singular.spoly{T}, 1},
        Z::Array{Singular.spoly{T}, 1}) where T
    @assert length(X) == length(Z)
    for i in 1:length(X)
        p = substitute_variable(p, X[i], Z[i])
    end
    return p
end

function Module{T <: Nemo.RingElem}(A::Array{Singular.spoly{T}, 2})
    R = Singular.parent(A[1])
    cols = [ Singular.svector(A[:, i]) for i in 1:size(A, 2) ]
    return Singular.Module(R, cols...)
end
