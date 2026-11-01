abstract type DistributionMatcher{X} end

match!(d::DistributionMatcher{X}, u::Vector{BIT}, dst::Vector{X}) where {X} = ABSTRACT(DM_TAG, "match!")
invert!(d::DistributionMatcher{X}, x::Vector{X}, u::Vector{BIT}) where {X} = ABSTRACT(DM_TAG, "invert!")
target_distribution(d::DistributionMatcher) = ABSTRACT(DM_TAG, "target_distribution")
info_rate(d::DistributionMatcher) = ABSTRACT(DM_TAG, "rate")
blocklength(d::DistributionMatcher) = ABSTRACT(DM_TAG, "blocklength")
messagelength(d::DistributionMatcher) = ABSTRACT(DM_TAG, "blocklength")

Base.show(io::IO, d::DistributionMatcher) = ABSTRACT(DM_TAG, "Base.show")
Base.parse(::Type{T}, s::AbstractString) where {T <: DistributionMatcher} = ABSTRACT(DM_TAG, T, "Base.parse")

abstract type PolarDistributionMatcher <: DistributionMatcher{BIT} end
include("polar_code/polar_code_distribution_matcher.jl")