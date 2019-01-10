import AbstractAlgebra
import Singular

#=
modified copy from Singular.jl (commit b277b7c) to make use_cache and
use_tensor_trick available
=#
function fres(id::Singular.sideal{T}, max_length::Int,
        method::String = "complete";
        use_cache::Bool = true, use_tensor_trick::Bool = false
        ) where T <: AbstractAlgebra.RingElem
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
    r, length = PrymGreen.id_fres(id.ptr, Cint(max_length), method, use_cache,
        use_tensor_trick, R.ptr)
    return Singular.sresolution{T}(R, Int(length), r)
end

#=
Apply the map x -> z to p where x is assumed to be a variable.
=#
function substitute_variable(p::Singular.spoly{T}, x::Singular.spoly{T},
        z::Singular.spoly{T}) where T <: AbstractAlgebra.RingElem
    R = Singular.parent(p)
    all(a -> Singular.parent(a) === R, [x, z]) || error("incompatible rings")
    index_var = findfirst(isequal(x), Singular.gens(R))
    index_var != nothing || error("x is not a variable")
    res_ptr = p_SubstPoly(p.ptr, index_var, z.ptr, R.ptr, R.ptr);
    return R(res_ptr)
end

#=
Successively apply the maps X[i] -> Z[i] to p where the X[i] are assumed to be
variables.
=#
function poly_substitute(p::Singular.spoly{T}, X::Array{Singular.spoly{T}, 1},
        Z::Array{Singular.spoly{T}, 1}) where T <: AbstractAlgebra.RingElem
    length(X) == length(Z) || error("incompatible lengths")
    for i in 1:length(X)
        p = substitute_variable(p, X[i], Z[i])
    end
    return p
end

function Module(A::Array{Singular.spoly{T}, 2}
        ) where T <: AbstractAlgebra.RingElem
    R = Singular.parent(A[1])
    all(a -> Singular.parent(a) === R, A) || error("incompatible rings")
    cols = [ Singular.vector(R, A[:, i]...) for i in 1:size(A, 2) ]
    return Singular.Module(R, cols...)
end

function Int(n::Singular.n_Z)
    R = Singular.parent(n)
    return n_Int(n.ptr, R.ptr)
end

function unwrap(n::AbstractAlgebra.RingElem)
    if typeof(n) <: Singular.n_unknown
        return Singular.libSingular.julia(n.ptr)
    end
    return n
end

Broadcast.broadcastable(n::AbstractAlgebra.RingElem) = Ref(n)
