using Distributions

#=
    Generate a vector of binary values (0 or 1) of length N, where each value is generated uniformly at random.
    The generated vector is stored in the provided destination vector dst.
    @param dst the destination vector to store the generated binary values
    @return the destination vector dst containing the generated binary values
=#
function gen_uniform_binary_vec!(dst::Vector{T}) where {T <: Union{Bool, BIT}}
    @inbounds for i in eachindex(dst)
        dst[i] = T(rand(Bool))
    end
    return dst
end

#=
    Generate a vector of binary values (0 or 1) of length N, where each value is generated uniformly at random.
    The generated vector is stored in the provided destination vector dst.
    @param N the length of the vector to generate
    @return the destination vector dst containing the generated binary values
=#
gen_uniform_binary_vec(::Type{T}, N::NAT) where {T <: Union{Bool, BIT}} = gen_uniform_binary_vec!(zeros(T, N))

#=
    Generate a vector of binary values (0 or 1) of length N, where each value is generated independently with probability p of being 1.
    The generated vector is stored in the provided destination vector dst.
    The function uses the Bernoulli distribution from the Distributions.jl package to generate the binary values.
    @param p the probability of generating a 1
    @param dst the destination vector to store the generated binary values
    @return the destination vector dst containing the generated binary values
=#
function gen_binary_vec!(p::REAL, dst::Vector{T}) where {T <: Union{Bool, BIT}}
    0 <= p <= 1 || ERROR(UTILS_VECT_GENERATOR_TAG, "$p is not a valid probability")
    p == 0.5 && return gen_uniform_binary_vec!(dst)
    sber = Bernoulli(p)
    @inbounds for i in eachindex(dst)
        dst[i] = T(rand(sber))
    end
    return dst
end

#=
    Generate a vector of binary values (0 or 1) of length N, where each value is generated independently with probability p of being 1.
    The function uses the Bernoulli distribution from the Distributions.jl package to generate the binary values.
    @param p the probability of generating a 1
    @param N the length of the vector to generate
    @return a vector containing the generated binary values
=#
gen_binary_vec(::Type{T}, N::NAT, p::REAL) where {T <: Union{Bool, BIT}} = gen_binary_vec!(p, zeros(T, N))

print_vector(v::Vector{T}) where {T} = begin
    s = ""
    for i in v
        s = s * "$i " 
    end
    return s
end 