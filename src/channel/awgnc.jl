using Distributions


struct AWGNC <: Channel{REAL, REAL}
    _sigma::REAL
    _underlying_dist::Normal{REAL}
    _last_gen_noise::Vector{REAL}
    _last_signal_pow::Vector{REAL}
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