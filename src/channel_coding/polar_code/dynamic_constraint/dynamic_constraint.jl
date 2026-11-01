abstract type DynamicConstraint end
(D::DynamicConstraint)(u::Vector{BIT}, select_index::Vector{Bool}, index::NAT) = ABSTRACT(PCDC_CONSTRAINT_TAG, D, "dynamic function")
known_at_decoding(D::DynamicConstraint, index::NAT) = ABSTRACT(PCDC_CONSTRAINT_TAG, "known_at_decoding")
Base.show(io::IO, D::DynamicConstraint) = ABSTRACT(PCDC_CONSTRAINT_TAG, D, "Base.show")
Base.parse(::Type{T}, s::AbstractString) where {T <: DynamicConstraint} = ABSTRACT(PCDC_CONSTRAINT_TAG, T, "Base.parse")
clear_cache(D::DynamicConstraint) = ABSTRACT(PCDC_CONSTRAINT_TAG, "clear_cache")

include("crc_constraint.jl")
include("row_merge_constraint.jl")
include("none_constraint.jl")
include("pc_constraint.jl")

#=
function gen_shaping_constraint(args)
    if args["shaping_constraint"] == "None"
        return NoneConstraint
    elseif args["shaping_constraint"] == "CRC"
        return get_crc_constraint(args["crc_width"], args["crc_poly"], ceil(NAT, args["shaping_rate"]*args["blocklength"]))
    end
end
function gen_frozen_constraint(args)
    ceil(NAT, rate*blocklength)
end
=#