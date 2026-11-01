none_contraint_tag = "None Constraint"

struct NoneConstraint <: DynamicConstraint end

(N::NoneConstraint)(u::Vector{BIT}, select_index::Vector{Bool}, index::NAT) = ERROR(none_contraint_tag, "This function should never be called (was, on index $index)")
known_at_decoding(D::NoneConstraint, index::NAT) = ERROR(none_contraint_tag, "This function should never be called (was, on index $index)")
Base.show(io::IO, D::NoneConstraint) = print(io, "None()")
Base.parse(::Type{T}, s::AbstractString) where {T <: NoneConstraint} = begin
    if s == "None()"
        return NoneConstraint
    else
        ERROR(none_contraint_tag, "input string cannot be parsed")
    end
end
clear_cache(D::NoneConstraint) = begin end
