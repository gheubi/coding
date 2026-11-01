using Distributions

#=
    An additive white Gaussian noise channel is a channel that adds Gaussian noise to the input signal. 
    The channel is defined by the following parameters:
        - sigma: the standard deviation of the noise
    The channel is implemented using the Normal distribution from the Distributions.jl package.
=#
struct BEC <: Channel{BIT, BIT}
    _erasure_prob::REAL                 # erasure probability of the channel
    _underlying_ber::Bernoulli{REAL}    # underlying Bernoulli distribution used to generate erasures
    _last_gen_noise::Vector{REAL}       # last generated noise power of the channel, which is the number of symbols that were flipped during the last call to transmit!()
    _last_signal_pow::Vector{REAL}      # last signal power of the channel, which is the number of symbols that were transmitted during the last call to transmit!()
    function BEC(p::REAL)
        0 <= p <= 1|| ERROR(BEC_TAG, "$p is not a valid erasure probability")
        new(p, Bernoulli(p), zeros(REAL, NDEFAULT), zeros(REAL, NDEFAULT))
    end
end

transmit!(W::BEC, src::AbstractVector{BIT}, dst::AbstractVector{BIT}) =  begin
    tid = default_tid()
    W._last_gen_noise[tid] = 0
    W._last_signal_pow[tid] = 0
    @inbounds for i in eachindex(src)
        W._last_signal_pow[tid] += src[i]
        dst[i] = rand(W._underlying_ber) ? ERASURE : src[i]
        if dst[i] == ERASURE
            W._last_gen_noise[i] += 1
        end
    end
    return dst
end
transmit_one_symbol(W::BEC, x::BIT) = begin
    rand(W._underlying_ber) ? ERASURE : x
end

(W::BEC)(y::BIT, x::BIT) = begin
    x == 0 || x == 1 || ERROR(BEC_TAG, "BEC(y|x) : x ($x) must be in {0, 1}")
    y == 0 || y == 1 || y == ERASURE || ERROR(BEC_TAG, "BEC(y|x) : y ($y) must be in {0, 1, ?}")
    return y == ERASURE ? W._erasure_prob : y == x ? (1 - W._erasure_prob) : 0
end
Base.show(io::IO, W::BEC) = print(io, "BEC(erasure_prob=$(W._erasure_prob))")
Base.parse(::Type{T}, s::AbstractString) where {T<:BEC} = BEC(parse(REAL, match(r"BEC\(erasure_prob=(.*)\)", s).captures[1]))


last_gen_noise_pow(W::BEC) = W._last_gen_noise[default_tid()]
last_signal_pow(W::BEC) = W._last_signal_pow[default_tid()]