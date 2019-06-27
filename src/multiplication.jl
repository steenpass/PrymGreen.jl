import OffsetArrays

function max_values_per_row(g::Int)
    g == 8 && return 11
    g == 10 && return 75
    (g % 2 == 0 && g > 10) || error("invalid input")
    return (g-5)*binomial(div(g, 2), 3) + (g-7)*binomial(div(g, 2), 2)
end

#=
Return the maximal p::Entry_t such that (p-1)+N*(p^2-p) <= 2^E-1.
Works for g in 8:2:50.
=#
function max_modulus(g::Int)
    E = sizeof(PrymGreen.Arith_t)*8
    N = max_values_per_row(g)
    return floor(Entry_t, 0.5+sqrt(2.0^E/N))
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

h_shift(h, v, f) = binomial(h+v+f-3, h-1)
v_shift(h, v, f) = binomial(h+v+f-3, v-1)

function horizontal_blocks(g::Int)
    hblocks = [ n for n in (div(g, 2)-2):(g-5) ]
    push!(hblocks, g-7)
    return hblocks
end

function prym_green_size(g::Int)
    # sum up v_shift's of all blocks in the last column
    res = 3*binomial(g-5, div(g, 2)-4)
    res += g*binomial(g-4, div(g, 2)-4)
    return res
end

function count_values_block(h::Int, v::Int, f::Int)
    (h < 1 || v < 1 || f < 1) && return 0
    return binomial(h+v+f-3, f-1)
end

function count_values_row(g::Int, f::Int)
    n_values = 0
    v = div(g, 2)-f-1
    hblocks = horizontal_blocks(g)
    h = 1
    while h < size(hblocks, 1)
        n_values += hblocks[h]*count_values_block(h, v, f)
        n_values += hblocks[h]*count_values_block(h, v-1, f+1)
        h += 1
    end
    # h == size(hblocks, 1)
    n_values += hblocks[h]*count_values_block(h, v, f)
    return n_values
end

MultData = OffsetArrays.OffsetArray{Array{Tuple{Msize_t, Msize_t, Int8}, 1}, 1}

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

# y = A*x
function write_multiply_koszul_block(g::Int, h::Int, v::Int, f::Int,
        p::Entry_t)
    out = "static void multiply_koszul_block_g" * string(g) * "_"
    out *= string(h) * "_" * string(v) * "_" * string(f)
    out *= "(arith_t *y, arith_t *A, arith_t *x)\n{\n"
    D = multiplication_data(h, v, f)
    for i in axes(D, 1)
        offset = Int128(count(x -> (x[3] < 0), D[i]) * (p^2-p))
        out2 = "    y[" * string(i) * "] += " * string(offset) * Arith_t_suffix
        for (j, p, s) in D[i]
            out2 *= (s > 0 ? "+" : "-") *
                    "A[" * string(p) * "]*x[" * string(j) * "]"
        end
        out *= out2 * ";\n"
    end
    out *= "}\n\n"
    return out
end

function call_multiply_koszul_block(g::Int, h::Int, v::Int, f::Int,
        offset_y::Int, offset_A::Int, offset_x::Int)
    out = "    multiply_koszul_block_g" * string(g) * "_"
    out *= string(h) * "_" * string(v) * "_" * string(f)
    out *= "(y+" * string(offset_y)
    out *= ", A+" * string(offset_A)
    out *= ", x+" * string(offset_x) * ");\n"
    return out
end

function write_multiply_koszul_row(g::Int, f::Int)
    out = "static void multiply_koszul_row_g" * string(g) * "_f" * string(f)
    out *= "(arith_t *y, arith_t *A, arith_t *x)\n{\n"
    v = div(g, 2)-f-1
    offsets = [ 0, 0, 0 ]   # offsets for y, A, x
    hblocks = horizontal_blocks(g)
    for h in 1:size(hblocks, 1)
        for j in 1:hblocks[h]
            out *= call_multiply_koszul_block(g, h, v, f, offsets...)
            offsets[2] += count_values_block(h, v, f)
            if h != size(hblocks, 1)
                offsets[1] += v_shift(h, v, f)
                out *= call_multiply_koszul_block(g, h, v-1, f+1, offsets...)
                offsets[2] += count_values_block(h, v-1, f+1)
                offsets[1] = 0
            end
            offsets[3] += h_shift(h, v, f)
        end
    end
    out *= "}\n\n"
    return out
end

function call_multiply_koszul_row(g::Int, f::Int, offset_y::Int, offset_A::Int)
    out = "    multiply_koszul_row_g" * string(g) * "_f" * string(f)
    out *= "(y+" * string(offset_y)
    out *= ", A+" * string(offset_A)
    out *= ", x);\n"
    return out
end

function write_multiply_matrix_loop(g::Int)
    out = "static void multiply_matrix_loop_g" * string(g)
    out *= "(arith_t *y, arith_t *A, arith_t *x)\n{\n"
    offsets = [ 0, 0 ]   # offsets for y, A
    shift_f2 = [ v_shift(div(g, 2)-1, div(g, 2)-3, 2), count_values_row(g, 2) ]
    shift_f3 = [ v_shift(div(g, 2)-1, div(g, 2)-4, 3), count_values_row(g, 3) ]
    for k in 1:3
        out *= call_multiply_koszul_row(g, 2, offsets...)
        offsets += shift_f2
    end
    for k in 1:g
        out *= call_multiply_koszul_row(g, 2, offsets...)
        offsets += shift_f2
        out *= call_multiply_koszul_row(g, 3, offsets...)
        offsets += shift_f3
    end
    out *= "}\n\n"
    return out
end

function write_recurrence_sequence(g::Int)
    out = "msize_t recurrence_sequence_g" * string(g)
    out *= "(arith_t **seq, arith_t *A, arith_t *v, msize_t index)\n{\n"
    pg_size = prym_green_size(g)
    size_v = pg_size*sizeof(Arith_t)
    out *= "    arith_t *v_a = (arith_t *)malloc(" * string(size_v) * ");\n"
    out *= "    arith_t *v_b = (arith_t *)malloc(" * string(size_v) * ");\n"
    out *= "    memcpy(v_b, v, " * string(size_v) * ");\n"
    length = 2*pg_size
    out *= "    // this memory block will be handed over to the calling "
    out *= "function:\n"
    out *= "    *seq = (arith_t *)malloc(" * string(length*sizeof(Arith_t))
    out *= ");\n"
    out *= "    (*seq)[0] = v[index];\n"
    out *= "    for (msize_t i = 1; i < " * string(length) * "; i++) {\n"
    out *= "        memset(v_a, 0, " * string(size_v) * ");\n"
    out *= "        multiply_matrix_loop_g" * string(g) * "(v_a, A, v_b);\n"
    out *= "        (*seq)[i] = v_a[index];\n"
    out *= "        memcpy(v_b, v_a, " * string(size_v) * ");\n"
    out *= "    }\n"
    out *= "    free(v_a);\n"
    out *= "    free(v_b);\n"
    out *= "    return " * string(length) * ";\n"
    out *= "}\n\n"
end

function file_for_multiplication_code(g::Int)
    dir = joinpath(realpath(joinpath(@__DIR__, "..", "src", "lib")), "tmp")
    mkpath(dir)
    file = joinpath(dir, "prym_green_g" * string(g) * ".c")
    return file
end

function write_multiplication_code(g::Int, p::Entry_t)
    file = file_for_multiplication_code(g)
    io = open(file, "w")
    write(io, "#include <string.h>\n")
    write(io, "#include \"../prym_green_types.h\"\n\n")
    B = needed_blocks(g)
    for (h, v, f) in B
        write(io, write_multiply_koszul_block(g, h, v, f, p))
    end
    write(io, write_multiply_koszul_row(g, 2))
    write(io, write_multiply_koszul_row(g, 3))
    write(io, write_multiply_matrix_loop(g))
    write(io, write_recurrence_sequence(g))
    close(io)
end

function compile_multiplication_code(g::Int)
    file = file_for_multiplication_code(g)
    o_file = file[1:(end-1)] * "o"
    so_file = file[1:(end-1)] * "so"
    run(`gcc -c -fPIC -Wall -o $o_file $file`)
    run(`gcc -shared -fPIC $o_file -o $so_file`)
end
