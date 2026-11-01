#=
    Compute the weigt of integer i represented using n bits
=#
function weight(i::NAT, n::NAT)
    acc = 0
    for j ∈ 1:n
        acc += i&1
        i >>>= 1
    end
    return acc
end
#=
    Compute the weight of integer i as the number of 
    ones in its binary representation
=#
function weight(i::NAT)
    acc = 0
    while i != 0
        acc += i&1
        i >>>= 1
    end
    return acc
end
#=
    Compute the weight of a boolean vector, i.e.
    the number of true values in the vector
=#
function weight(v::Vector{Bool})
    acc = 0
    for b ∈ v
        if b 
            acc += 1
        end
    end
    return acc
end
#=
    Compute the weight of a vector of integers, i.e.
    the number of non-zero values in the vector
=#
function weight(v::Vector{<:Integer})
    acc = 0
    for b ∈ v
        acc += b != 0
    end
    return acc
end