import OffsetArrays

function max_values_per_row(g::Int)
    g == 8 && return 11
    g == 10 && return 75
    (g % 2 == 0 && g > 10) || error("invalid input")
    return (g-5)*binomial(div(g, 2), 3) + (g-7)*binomial(div(g, 2), 2)
end

function max_modulus(g::Int)
    E = sizeof(PrymGreen.Arith_t)*8
    N = max_values_per_row(g)
    return ceil(Entry_t, sqrt(2.0^E/N))
end

function needed_blocks(g::Int)
    (g % 2 == 0 && g > 6) || error("invalid input")
    (v, f) = (div(g, 2)-3, 2)
    B1 = [ (h, v, f) for h in 1:(div(g, 2)-1) ]
    (v, f) = (div(g, 2)-4, 3)
    B2 = [ (h, v, f) for h in 1:(div(g, 2)-1) if v > 0 ]
    (v, f) = (div(g, 2)-5, 4)
    B3 = [ (h, v, f) for h in 1:(div(g, 2)-2) if v > 0 ]   # !
    return vcat(B1, B2, B3)
end

MultData = OffsetArrays.OffsetArray{Array{Tuple{Msize_t, Msize_t, Int8}, 1}, 1}

h_shift(h, v, f) = binomial(h+v+f-3, h-1)
v_shift(h, v, f) = binomial(h+v+f-3, v-1)

function f_multiply(D::MultData, h::Int, v::Int, s::Int,
        i::Msize_t, j::Msize_t, p::Msize_t)
    size = h_shift(h, v, 1)   # == v_shift(h, v, 1)
    for k in 0:(size-1)
        push!(D[i+k], (j+k, p, s))
    end
end

function h_multiply(D::MultData, v::Int, f::Int, s::Int,
        i::Msize_t, j::Msize_t, p::Msize_t)
    size = v_shift(1, v, f)
    for k in 0:(size-1)
        push!(D[i+k], (j, p+size-k-1, s))
        s *= -1
    end
end

function v_multiply(D::MultData, h::Int, f::Int, s::Int,
        i::Msize_t, j::Msize_t, p::Msize_t)
    size = h_shift(h, 1, f)
    for k in 0:(size-1)
        push!(D[i], (j+k, p+k, s))
    end
end

function multiplication_data_recursion(D::MultData, h::Int, v::Int, f::Int,
        s::Int, i::Msize_t, j::Msize_t, p::Msize_t)
    (h < 1 || v < 1 || f < 1) && return
    (f == 1)           && return f_multiply(D, h, v, s, i, j, p)
    (h == 1 && f == 2) && return h_multiply(D, v, f, s, i, j, p)
    (v == 1)           && return v_multiply(D, h, f, s, i, j, p)
    multiplication_data_recursion(D, h-1, v, f, s, i, j, p)
    j += Msize_t(h_shift(h-1, v, f))
    f_p = p+Msize_t(binomial(h+v+f-4, f-1))
    multiplication_data_recursion(D, h, v, f-1, s, i, j, f_p)
    i += Msize_t(v_shift(h, v, f-1))
    s *= (f%2)*2-1
    multiplication_data_recursion(D, h, v-1, f, s, i, j, p)
end

function multiplication_data(h::Int, v::Int, f::Int)
    A = [ eltype(MultData)(undef, 0) for i in 1:v_shift(h, v, f) ]
    D = OffsetArrays.OffsetArray(A, -1)
    n = binomial(h+f-2, f-1)
    for x in D
        sizehint!(x, n)
    end
    (i, j, p) = (Msize_t(0), Msize_t(0), Msize_t(0))
    multiplication_data_recursion(D, h, v, f, 1, i, j, p)
    for x in D
        sort!(x)
    end
    return D
end

