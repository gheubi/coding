bit_rev_op(i::NAT, n::NAT) = begin
    rev = 0
    ire = i-1
    for _ in 1:n
        rev <<= 1
        rev |= ire&1
        ire >>= 1
    end
    rev + 1
end


function polar_transform!(u::Vector{BIT})
    N = length(u)
    @check ispow2(N) || ERROR(POLAR_TRANSFORM, "$N should be a power of two")
    n = convert(NAT, log2(N))

    @inbounds for i in 1:N
        j = bit_rev_op(i, n)
        if j > i
            u[i], u[j] = u[j], u[i]
        end
    end
    @inbounds for i in 1:n
        steps = 1 << (i-1)
        chunk = N ÷ steps
        halfchunk = chunk >> 1
        for j in 1:steps
            floor = chunk*(j-1)
            for k in 1:halfchunk
                u[floor + k] ⊻= u[floor + halfchunk + k]
                u[floor + k] %= 2
            end
        end
    end
    return u
end