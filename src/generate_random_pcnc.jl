import AbstractAlgebra
import Nemo
import Singular

#=
Check that n is the multiplicative order of z, given that z^n = 1.
=#
function has_multiplicative_order(z::AbstractAlgebra.RingElem, n::Int)
    if z == 1 && n > 1
        return false
    end
    fac = Nemo.factor(Nemo.FlintZZ(n))
    for (p, e) in fac
        if z^div(n, Int(p)) == 1
            return false
        end
    end
    return true
end

#=
Choose a uniformly random prime p from the interval [a, b] such that there
exists a primitive n-th root of unity in Z/pZ, that is, such that n divides
p-1.
=#
function random_prime_with_primitive_root_of_unity(a::Int, b::Int, n::Int,
        rng::AbstractRNG)
    from = div(a+n-2, n)   # = ceil((a-1)/n)
    to = div(b-1, n)   # = floor((b-1)/n)
    interval = from:to
    contains_prime = false
    for k in interval
        p = n*k+1
        if Nemo.isprime(Nemo.fmpz(p))
            contains_prime = true
            break
        end
    end
    contains_prime || error("interval to small")
    k = rand(rng, interval)
    p = n*k+1
    while !Nemo.isprime(Nemo.fmpz(p))
        k = rand(rng, interval)
        p = n*k+1
    end
    p
end

#=
Choose a uniformly random prime p from the interval [a, b] such that there
exists a primitive n-th root of unity in R := Z/pZ and choose a uniformly
random n-th root of unity z in R. Return (R, z).
=#
function random_primitive_root_of_unity(a::Int, b::Int, n::Int,
        rng::AbstractRNG)
    n > 0 || error("n must be positive")
    p = random_prime_with_primitive_root_of_unity(a, b, n, rng)
    k = div(p-1, n)
    # R = Nemo.ResidueRing(Nemo.FlintZZ, p)
    R = Singular.Fp(p)
    z = R(rand(rng, 1:(p-1)))^k
    while !has_multiplicative_order(z, n)
        z = R(rand(rng, 1:(p-1)))^k
    end
    (R, z)
end

function are_distinct_points_of_P1{T <: AbstractAlgebra.RingElem}(
        M::Array{T, 2})
    @assert size(M, 1) == 2
    for i in 2:size(M, 2)
        for j in 1:(i-1)
            if M[1, j]*M[2, i]-M[1, i]*M[2, j] == 0
                return false
            end
        end
    end
    return true
end

function random_field_elements(rng::AbstractRNG, R::Nemo.NmodRing, dims::Dims)
    @assert Nemo.isprime(R.n)
    return R.(rand(rng, 1:R.n, dims))
end

function random_field_elements(rng::AbstractRNG, R::Singular.N_ZpField,
        dims::Dims)
    return R.(rand(rng, 1:Int(Singular.characteristic(R)), dims))
end

function random_distinct_points_of_P1(R::AbstractAlgebra.Ring, n::Int,
        rng::AbstractRNG)
    M = random_field_elements(rng, R, (2, n))
    while !are_distinct_points_of_P1(M)
        M = random_field_elements(rng, R, (2, n))
    end
    return M
end

function canonical_multipliers{T <: AbstractAlgebra.RingElem}(P::Array{T, 2},
        Q::Array{T, 2})
    @assert size(P, 1) == 2
    @assert size(Q, 1) == 2
    g = size(P, 2)
    R = AbstractAlgebra.parent(P[1])
    S, (x_0, x_1) = Singular.PolynomialRing(R, ["x_0", "x_1"])
    zeros_P = [ (P[1, i]*x_1-P[2, i]*x_0) for i in 1:g ]
    zeros_Q = [ (Q[1, i]*x_1-Q[2, i]*x_0) for i in 1:g ]
    quadrics = [ zeros_P[i]*zeros_Q[i] for i in 1:g ]
    sections = [ prod(quadrics[1:g .!= i]) for i in 1:g ]
    AP = [ poly_substitute(sections[i], [x_0, x_1], S.(P[:, i])) for i = 1:g ]'
    AQ = [ poly_substitute(sections[i], [x_0, x_1], S.(Q[:, i])) for i = 1:g ]'
    A = Singular.coeff.(vcat(AP, AQ), 0)
    return unwrap.(A)
end

function change_multiplier{T <: AbstractAlgebra.RingElem}(A::Array{T, 2}, r::T)
    @assert size(A, 1) == 2
    A[2, :] .*= r
    return A
end

function linear_series_from_multipliers{T <: AbstractAlgebra.RingElem}(
        P::Array{T, 2}, Q::Array{T, 2}, A::Array{T, 2})
    R = AbstractAlgebra.parent(P[1])
    S, X = Singular.PolynomialRing(R, ["x_0", "x_1"])
    g = size(P, 2)
    B = Singular.MaximalIdeal(S, 2*g-2)
    d = Singular.ngens(B)
    MP = [ A[2, i]*poly_substitute(B[j], X, S.(P[:, i])) for i = 1:g, j = 1:d ]
    MQ = [ A[1, i]*poly_substitute(B[j], X, S.(Q[:, i])) for i = 1:g, j = 1:d ]
    sy = Singular.syz(Module(MP-MQ))
    l = Singular.ngens(sy)
    gens = [ sum([ B[j] for j = 1:d ] .* Array(sy[i])) for i = 1:l ]
    map = Singular.Ideal(S, gens)
    return map
end

#=
Compute a random Prym canonical nodal curve of genus g and level l.
=#
function random_PCNC(g::Int, l::Int, rng::AbstractRNG)
    # 268435399 is the limit for Singular.Fp()
    (R, r) = random_primitive_root_of_unity(2, 268435399, l, rng)
    P = random_distinct_points_of_P1(R, g, rng)
    Q = random_distinct_points_of_P1(R, g, rng)
    A = canonical_multipliers(P, Q)
    A = change_multiplier(A, r)
    s = linear_series_from_multipliers(P, Q, A)
    n = Singular.ngens(s)
    @assert n == g-1
    T, = Singular.PolynomialRing(R, [ "t_$i" for i = 0:(n-1) ])
    K = Singular.kernel(T, s)
    return K
end
