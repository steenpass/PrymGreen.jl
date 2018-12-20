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
    r, length = PrymGreen.id_fres(id.ptr, Cint(max_length), method, R.ptr,
            use_cache, use_tensor_trick)
    return Singular.sresolution{T}(R, Int(length), r)
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

function Module{T <: AbstractAlgebra.RingElem}(A::Array{Singular.spoly{T}, 2})
    R = Singular.parent(A[1])
    cols = [ Singular.vector(R, A[:, i]...) for i in 1:size(A, 2) ]
    return Singular.Module(R, cols...)
end

function n_Int(n::Singular.libSingular.number, R::Singular.libSingular.coeffs)
    icxx"""n_Int($n, $R);"""
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
