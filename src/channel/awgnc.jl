using Distributions

#=
    An additive white Gaussian noise channel is a channel that adds Gaussian noise to the input signal. 
    The channel is defined by the following parameters:
        - sigma: the standard deviation of the noise
    The channel is implemented using the Normal distribution from the Distributions.jl package.
=#
struct AWGNC <: Channel{REAL, REAL}
    _sigma::REAL                        # standard deviation of the noise
    _underlying_dist::Normal{REAL}      # underlying normal distribution used to generate noise
    _last_gen_noise::Vector{REAL}       # last generated noise power of the channel, which is the number of symbols that were flipped during the last call to transmit!()
    _last_signal_pow::Vector{REAL}      # last signal power of the channel, which is the number of symbols that were transmitted during the last call to transmit!()
    function AWGNC(sigma::REAL)
        0 <= sigma || ERROR(AWGNC_TAG, "$sigma is not a valid standard deviation")
        new(sigma, Normal(0.0, sigma), zeros(NDEFAULT), zeros(NDEFAULT))
    end
end

snr_db_to_sigma(snr::REAL) = sqrt(1/value_from_db(snr))

transmit!(W::AWGNC, src::AbstractVector{REAL}, dst::AbstractVector{REAL}) =  begin
    n::REAL = 0.0
    tid = default_tid()
    W._last_gen_noise[tid] = 0.0
    W._last_signal_pow[tid] = 0.0
    @inbounds for i in eachindex(src)
        n = rand(W._underlying_dist)
        dst[i] = src[i] + n
        W._last_gen_noise[tid] += n^2
        W._last_signal_pow[tid] += src[i]^2
    end
    return dst
end
transmit_one_symbol(W::AWGNC, x::REAL) = begin
    x + rand(W._underlying_dist)
end
(W::AWGNC)(y::REAL, x::REAL) = begin
    return exp(-(y-x)^2/2(W._sigma^2))/(sqrt(2*pi*(W._sigma^2)))  
end
Base.show(io::IO, W::AWGNC) = print(io, "AWGNC(sigma=$(W._sigma))")
Base.parse(::Type{T}, s::AbstractString) where {T<:AWGNC} = AWGNC(parse(REAL, match(r"AWGNC\(sigma=(.*)\)", s).captures[1]))

last_gen_noise_pow(W::AWGNC) = W._last_gen_noise[default_tid()]
last_signal_pow(W::AWGNC) = W._last_signal_pow[default_tid()]