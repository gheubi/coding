
abstract type Code{U,X,Y} end

encode!(c::Code{U,X,Y}, u::Vector{U}, dst::Vector{X}) where {U, X, Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "encode!")
decode!(c::Code{U,X,Y}, y::Vector{Y}, dst::Vector{U}) where {U, X, Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "decode!")
blocklength(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "blocklength")
messagelength(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "messagelength")
set_channel(c::Code{U,X,Y}, W::Channel{X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "set_channel")
current_channel(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "current_channel")
is_valid_codeword(c::Code{U,X,Y}, x::Vector{X}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "is_valid_codeword")
rate(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "rate")
shaping_rate(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "shaping_rate")
target_dist(c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(CHANNEL_CODE_TAG, c, "target_dist")

Base.show(io::IO, c::Code{U,X,Y}) where {U,X,Y} = ABSTRACT(channel_tag, c, "Base.show")
Base.parse(::Type{T}, s::AbstractString) where {T<:Code} = ABSTRACT(channel_tag, T, "Base.parse")

abstract type PolarCode{Y} <: Code{BIT, BIT, Y} end
include("polar_code/polar_code.jl")

