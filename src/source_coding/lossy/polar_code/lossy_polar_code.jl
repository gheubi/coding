#=
    A concrete implementation of a lossy polar code.
    This struct wraps around a PolarCode and provides methods for encoding, decoding, and calculating distortion.
    It also provides methods to retrieve the input distribution, rate, blocklength, and messagelength of the code.
    The encode! method encodes a vector of input bits into a vector of output bits using the underlying PolarCode's decode! method.
    The decode! method decodes a vector of output bits into a vector of reconstructed bits using the underlying PolarCode's encode! method.
    The distortion method calculates the Hamming distance between the original input bits and the reconstructed bits.
    The input_distribution method returns the target distribution of the underlying PolarCode.
    The rate, blocklength, and messagelength methods return the corresponding properties of the underlying PolarCode.
=#
struct ConcreteLossyPolarCode <: LossyPolarCode
    _pc::PolarCode
end

encode!(c::ConcreteLossyPolarCode, u::Vector{BIT}, dst::Vector{BIT}) = begin
    decode!(c._pc, u, dst)
end
decode!(c::ConcreteLossyPolarCode, x::Vector{BIT}, dst::Vector{BIT}) = begin
    encode!(c._pc, x, dst)
end
distortion(c::ConcreteLossyPolarCode, u::Vector{BIT}, û::Vector{BIT}) = begin
    acc = 0
    @inbounds for i in eachindex(u)
        if u[i] != û[i]
            acc += 1
        end
    end
    return acc
end
input_distribution(c::ConcreteLossyPolarCode) = target_dist(c._pc)
rate(c::ConcreteLossyPolarCode) = rate(c._pc)
blocklength(c::ConcreteLossyPolarCode) = blocklength(c._pc)
messagelength(c::ConcreteLossyPolarCode) = messagelength(c._pc)

get_lossy_polar_code(dic) = ConcreteLossyPolarCode(get_polar_code(dic))

Base.show(io::IO, c::ConcreteLossyPolarCode) = print(io, "$(c._pc)")
Base.parse(::Type{T}, s::AbstractString) where {T <: ConcreteLossyPolarCode} = ABSTRACT(LOSSY_POLAR_CODE_TAG, T, "Base.parse")

