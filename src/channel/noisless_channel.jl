struct NoiselessChannel{U} <: Channel{U,U} 
    last_signal_pow::Vector{REAL}
end

NoiselessChannel{U}() where {U} = NoiselessChannel{U}(zeros(REAL,NDEFAULT))
    

(W::NoiselessChannel{U})(y::U, x::U) where {U} = begin 
    y == x ? 1.0 : 0.0
end

transmit!(W::NoiselessChannel{U}, src::AbstractVector{U}, dst::AbstractVector{U}) where {U} = begin
    tid = default_tid()
    W.last_signal_pow[tid] = 0.0
    @inbounds for i in eachindex(src)
        W.last_signal_pow[tid] += src[i]^2
        dst[i] = src[i]
    end
    return dst
end

transmit_one_symbol(W::NoiselessChannel{U}, x::U) where {U} = x

last_signal_pow(W::NoiselessChannel{U}) where {U} = W.last_signal_pow[default_tid()]
last_gen_noise_pow(W::NoiselessChannel) = 0.0

Base.show(io::IO, W::NoiselessChannel{U}) where {U} = print(io, "NoiselessChannel")
Base.parse(::Type{T}, s::AbstractString) where {U,T <: NoiselessChannel{U}} = begin
    if s == "NoiselessChannel"
        return NoiselessChannel{U}()
    else
        ERROR(NOISELESS_CHANEL_TAG, "cannot parse input string")
    end
end

LLR(W::NoiselessChannel{BIT}, y::BIT) = 0.0