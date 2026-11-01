#=
    Compute the bit-reversal of an integer i with respect to n bits.
    The function takes an integer i and the number of bits n as input, and returns 1 + the integer obtained by reversing the order of the bits in the binary representation of i-1, considering only the least
    significant n bits.
    @warning : If i is greater than or equal to 2^n, the function will only consider the least significant n bits of i for the bit-reversal operation.
    @warning : The function computes the bit-reversal for i in [1, 2^n]. It does so by shifting i to [0, 2^n-1], reversing the bits, and then shifting back to [1, 2^n].
    @example : If i = 7, i-1 = 6 (binary 110) and n = 3, the bit-reversal of i is 3 (binary 011) + 1 = 4.
    @param i the integer to be bit-reversed
    @param n the number of bits to consider for the bit-reversal
    @return the bit-reversed integer
=#
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

#=
    Compute the polar transform of a vector of bits u in-place.
    The function takes a vector of bits u as input and applies the polar transform to it in-place. 
    The polar transform is a linear transformation that combines the bits in u according to the polar coding scheme.

    polar_transform(u) = uB_nF^{\otimes n} where B_n is the bit-reversal permutation matrix and F^{\otimes n} is the n-th Kronecker power of the basic polar transform matrix F = [1 0; 1 1].

    @warning : The length of the input vector u must be a power of two. If it is not, the function will throw an error.
    @warning : The function modifies the input vector u in-place, so the original values in u will be lost after the transformation.
    @warning : In the litterature, the polar transform is sometiems defined without the bit-reversal permutation. In this implementation, the polar transform includes the bit-reversal permutation.

    @param u the vector of bits to be transformed
    @return the transformed vector of bits u after applying the polar transform
=#
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