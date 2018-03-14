#=
Compute the multiplicative order of z 
=#
function multiplicative_order(z::Nemo.nmod, n::Int)
    a = z
    ord = 1
    while a != 1
        a = a*z
        ord = ord+1
    end
    ord
end

#=
Choose a uniformly random prime p from the interval [a, b] such that there
exists a primitive n-th root of unity in Z/pZ and choose a uniformly random
n-th root of unity z in Z/pZ.
=#
function random_primitive_root_of_unity(a::Entry_t, b::Entry_t, n::Int,
        rng::AbstractRNG)
    n > 0 || error("n must be positive")
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
    R = Nemo.ResidueRing(Nemo.FlintZZ, p)
    z = R(rand(rng, 1:(p-1)))^k
    while multiplicative_order(z, n) != n
        z = R(rand(rng, 1:(p-1)))^k
    end
    (Entry_t(p), Entry_t(z.data))
end
