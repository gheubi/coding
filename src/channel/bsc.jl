using Distributions

#=
    An binary symmetric channel is a channel that flips the input bits with a certain probability. 
    The channel is defined by the following parameters:
        - p: the flip probability of the channel
    The channel is implemented using the Bernoulli distribution from the Distributions.jl package.
=#
struct BSC <: Channel{BIT, BIT}
    _flip_prob::REAL                        # flip probability of the channel
    _underlying_ber::Bernoulli{REAL}        # underlying Bernoulli distribution used to generate flips
    _last_gen_noise::Vector{REAL}           # last generated noise power of the channel, which is the number of symbols that were flipped during the last call to transmit!()
    _last_signal_pow::Vector{REAL}          # last signal power of the channel, which is the number of symbols that were transmitted during the last call to transmit!()
    function BSC(p::REAL)
        0 <= p <= 0.5 || ERROR(BSC_TAG, "$p should be within [0, 0.5]")
        new(p, Bernoulli(p), zeros(REAL, NDEFAULT), zeros(REAL, NDEFAULT))
    end
end

transmit!(W::BSC, src::AbstractVector{BIT}, dst::AbstractVector{BIT}) =  begin
    tid = default_tid()
    W._last_gen_noise[tid] = 0
    W._last_signal_pow[tid] = 0
    @inbounds for i in eachindex(src)
        W._last_signal_pow[tid] += src[i]
        dst[i] = rand(W._underlying_ber) ? (src[i]⊻1) : src[i]
        if dst[i] != src[i]
            W._last_gen_noise[tid] += 1
        end 
    end
    return dst
end
transmit_one_symbol(W::BSC, x::BIT) = begin
    rand(W._underlying_ber) ? x⊻1 : x
end

(W::BSC)(y::BIT, x::BIT) = begin
    x == 0 || x == 1 || ERROR(BSC_TAG, "BSC(y|x) : x ($x) must be in {0, 1}")
    y == 0 || y == 1 || ERROR(BSC_TAG, "BSC(y|x) : y ($y) must be in {0, 1}")
    return y != x ? W._flip_prob : (1 - W._flip_prob)
end
Base.show(io::IO, W::BSC) = print(io, "BSC(flip_prob=$(W._flip_prob))")
Base.parse(::Type{T}, s::AbstractString) where {T<:BSC} = BSC(parse(REAL, match(r"BSC\(flip_prob=(.*)\)", s).captures[1]))

last_gen_noise_pow(W::BSC) = W._last_gen_noise[default_tid()]
last_signal_pow(W::BSC) = W._last_signal_pow[default_tid()]