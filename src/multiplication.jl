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

