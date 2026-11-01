#=
    A binary constellation channel is a channel that maps a binary input to a real output using a given constellation. 
    The channel is defined by the following parameters:
        - N: the number of symbols to transmit
        - sigma: the standard deviation of the noise
        - constellation: a tuple of two real numbers representing the constellation points for 0 and 1 respectively.
    The channel is implemented as a wrapper around an AWGN channel, where the binary input is mapped to the corresponding constellation point and then transmitted through the AWGN channel.
=#  
struct BinaryConstellation <: Channel{BIT, REAL}
    _w::AWGNC
    _buffers::Vector{Vector{REAL}}
    _const::NTuple{2, REAL}
    function BinaryConstellation(N::NAT, sigma::REAL, constellation::NTuple{2, REAL})
        buff = Vector{Vector{REAL}}(undef, NDEFAULT)
        @inbounds for i in 1:NDEFAULT
            buff[i] = zeros(REAL, N)
        end
        new(AWGNC(sigma), buff, constellation)
    end
end

(W::BinaryConstellation)(y::REAL, u::BIT) = W._w(y, W._const[u+1])

transmit!(W::BinaryConstellation, src::AbstractVector{BIT}, dst::AbstractVector{REAL}) = begin
    length(src) == length(dst) == length(W._buffers[1]) || ERROR(BIN_CONST_AWGNC_TAG, "dimension mismatch")

    v = W._buffers[default_tid()]
    @inbounds for i in eachindex(v)
        v[i] = W._const[src[i]+1]
    end
    transmit!(W._w, v, dst)
end

transmit_one_symbol(W::BinaryConstellation, x::BIT) = begin
    transmit_one_symbol(W._w, W._const[x+1])
end

#=
    Show the channel information. 
    For OOK and BPSK constellations, the channel is shown as OOK or BPSK respectively. 
    For other constellations, the channel is shown as BinaryCons.
=#
Base.show(io::IO, W::BinaryConstellation) = begin 
    if W._const[1] == 0
        print(io, "OOK(sigma=$(W._w._sigma);p=$(1/(W._const[2]^2));N=$(length(W._buffers[1])))")
    elseif W._const[1] == 1 && W._const[2] == -1
        print(io, "BPSK(sigma=$(W._w._sigma);N=$(length(W._buffers[1])))")
    else
        print(io, "BinaryCons(sigma=$(W._w._sigma);constellation=[$(W._const[1]), $(W._const[2])];N=$(length(W._buffers[1])))")
    end
end

Base.parse(::Type{T}, s::AbstractString) where {T<:BinaryConstellation} = begin 
    if startswith("OOK(", s) 
        dic = parse_to_dic(s, "OOK", ["sigma", "p", "N"])
        return OOK(parse(NAT, dic["N"]), parse(REAL, dic["sigma"]), parse(REAL, dic["p"]))
    elseif startswith("BPSK(", s)
        dic = parse_to_dic(s, "BPSK", ["sigma", "N"])
        return BPSK(parse(NAT, dic["N"]), parse(REAL, dic["sigma"]))
    elseif startswith("BinaryCons(")
        dic = parse_to_dic(s, "BinaryCons", ["sigma", "constellation", "N"])
        con = NTuple{2, REAL}(parse_vector(REAL, dic["constellation"]))
        return BinaryConstellation(parse(NAT, dic["N"]), parse(NAT, dic["sigma"]), con)
    else
        ERROR(BIN_CONST_AWGNC_TAG, "Unrecongised string $s")
    end
end

last_gen_noise_pow(W::BinaryConstellation) = last_gen_noise_pow(W._w)
last_signal_pow(W::BinaryConstellation) = last_signal_pow(W._w)

#=
    Create a BPSK constellation channel with N symbols and noise standard deviation sigma.
=#
BPSK(N::NAT, sigma::REAL) = BinaryConstellation(N, sigma, (1.0, -1.0))
#=
    Create an OOK constellation channel with N symbols, noise standard deviation sigma and source probability source_p. 
    The constellation points are (0, sqrt(1/source_p)) to ensure that the average power of the transmitted symbols is 1.
=#
OOK(N::NAT, sigma::REAL, source_p::REAL) = begin 
    0<= source_p <= 1 || ERROR(BIN_CONST_AWGNC_TAG, "$source_p is not a valid probabiltiy for a OOK constellation") 
    BinaryConstellation(N, sigma, (0.0, sqrt(1/source_p)))
end

