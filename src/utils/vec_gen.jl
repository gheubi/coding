using Distributions

function gen_uniform_binary_vec!(dst::Vector{T}) where {T <: Union{Bool, BIT}}
    @inbounds for i in eachindex(dst)
        dst[i] = T(rand(Bool))
    end
    return dst
end

gen_uniform_binary_vec(::Type{T}, N::NAT) where {T <: Union{Bool, BIT}} = gen_uniform_binary_vec!(zeros(T, N))

function gen_binary_vec!(p::REAL, dst::Vector{T}) where {T <: Union{Bool, BIT}}
    0 <= p <= 1 || ERROR(UTILS_VECT_GENERATOR_TAG, "$p is not a valid probability")
    p == 0.5 && return gen_uniform_binary_vec!(dst)
    sber = Bernoulli(p)
    @inbounds for i in eachindex(dst)
        dst[i] = T(rand(sber))
    end
    return dst
end

gen_binary_vec(::Type{T}, N::NAT, p::REAL) where {T <: Union{Bool, BIT}} = gen_binary_vec!(p, zeros(T, N))

print_vector(v::Vector{T}) where {T} = begin
    s = ""
    for i in v
        s = s * "$i " 
    end
    return s
end 