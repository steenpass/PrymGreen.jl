#=
Check that n is the multiplicative order of z, given that z^n = 1.
=#
function has_multiplicative_order(z::Nemo.nmod, n::Int)
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
    R = Nemo.ResidueRing(Nemo.FlintZZ, p)
    z = R(rand(rng, 1:(p-1)))^k
    while !has_multiplicative_order(z, n)
        z = R(rand(rng, 1:(p-1)))^k
    end
    (R, Entry_t(z.data))
end

function are_distinct_points_of_P1(M::Array{Nemo.nmod, 2})
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

function random_distinct_points_of_P1(R::Nemo.NmodRing, n::Int,
        rng::AbstractRNG)
    M = R.(rand(rng, 1:R.n, (2, n)))
    while !are_distinct_points_of_P1(M)
        M = R.(rand(rng, 1:R.n, (2, n)))
    end
    M
end

function canonical_multipliers(P::Array{Nemo.nmod, 2}, Q::Array{Nemo.nmod, 2})
    @assert size(P, 1) == 2
    @assert size(Q, 1) == 2
    g = size(P, 2)
    R = Nemo.parent(P[1])
    S, (x_0, x_1) = Singular.PolynomialRing(R, ["x_0", "x_1"])
    zeros_P = [ (P[1, i]*x_1-P[2, i]*x_0) for i in 1:g ]
    zeros_Q = [ (Q[1, i]*x_1-Q[2, i]*x_0) for i in 1:g ]
    quadrics = [ zeros_P[i]*zeros_Q[i] for i in 1:g ]
    sections = [ prod(quadrics[1:g .!= i]) for i in 1:g ]
    AP = [ poly_substitute(sections[i], [x_0, x_1], S.(P[:, i])) for i = 1:g ]'
    AQ = [ poly_substitute(sections[i], [x_0, x_1], S.(Q[:, i])) for i = 1:g ]'
    A = Singular.coeff.(vcat(AP, AQ), 0)
    return (x -> Singular.libSingular.julia(x.ptr)).(A)
end

#=
Compute a random Prym canonical nodal curve of genus g and level l.
=#
function random_PCNC(g::Int, l::Int, rng::AbstractRNG)
    (R, z) = random_primitive_root_of_unity(2, 2147483647, l, rng)
    P = random_distinct_points_of_P1(R, g, rng)
    Q = random_distinct_points_of_P1(R, g, rng)
    A = canonical_multipliers(P, Q)
end
