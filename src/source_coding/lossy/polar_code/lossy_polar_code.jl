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

