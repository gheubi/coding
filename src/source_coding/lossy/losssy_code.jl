abstract type LossyCode{U,X,Y} end

encode!(c::LossyCode{U,X,Y}, u::Vector{U}, dst::Vector{X}) where {U,X,Y} = ABSTRACT(LOSSY_CODE_TAG, "encode!")
decode!(c::LossyCode{U,X,Y}, x::Vector{X}, dst::Vector{Y}) where {U,X,Y} = ABSTRACT(LOSSY_CODE_TAG, "decode!")
distortion(c::LossyCode{U,X,Y}, u::Vector{U}, û::Vector{Y}) where {U,X,Y} = ABSTRACT(LOSSY_CODE_TAG, "distortion")
input_distribution(c::LossyCode) = ABSTRACT(LOSSY_CODE_TAG, "input_distribution")
rate(c::LossyCode) = ABSTRACT(LOSSY_CODE_TAG, "rate")
blocklength(c::LossyCode) = ABSTRACT(LOSSY_CODE_TAG, "blocklength")
messagelength(c::LossyCode) = ABSTRACT(LOSSY_CODE_TAG, "messagelength")

Base.show(io::IO, c::LossyCode) = ABSTRACT(LOSSY_CODE_TAG, c, "Base.show")
Base.parse(::Type{T}, s::AbstractString) where {T <: LossyCode} = ABSTRACT(LOSSY_CODE_TAG, T, "Base.parse")

abstract type LossyPolarCode <: LossyCode{BIT,BIT,BIT} end
include("polar_code/lossy_polar_code.jl")