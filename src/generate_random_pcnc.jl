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
