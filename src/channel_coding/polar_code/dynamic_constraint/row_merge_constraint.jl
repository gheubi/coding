struct RowMergeConstraint <: DynamicConstraint
    _n_merges::NAT
    _merges::Vector{NTuple{2, NAT}}
    _previous_value_for::Dict{NAT, NAT}
end

function get_row_merge_constraint(merges::Vector{NTuple{2, NAT}})
    n_merges = length(merges)
    previous_value_for = Dict{NAT, NAT}()
    for t in merges
        t[2] >= t[1] || ERROR(PCDC_CONSTRAINT_ROWMERGE_TAG, "info bit must be after a dynamically frozen bit")
        previous_value_for[t[2]] = t[1]
    end
    return RowMergeConstraint(n_merges, merges, previous_value_for)
end

(R::RowMergeConstraint)(u::Vector{BIT}, select_index::Vector{Bool}, index::NAT) = begin
    return u[R._previous_value_for[index]]
end

Base.show(io::IO, R::RowMergeConstraint) = print(io, "RowMerge(n_merges=$(R._n_merges);merges=$(R._merges))")

Base.parse(::Type{T}, s::AbstractString) where {T <: RowMergeConstraint} = begin
    startswith(s, "RowMerge(") || ERROR(PCDC_CONSTRAINT_ROWMERGE_TAG, "input string cannot be parsed")

    dic = parse_to_dic(s, "RowMerge", ["n_merges", "merges"])
    n_merges = parse(NAT, dic["n_merges"])

    bad_merges = parse.(NAT, split(replace(replace(replace(replace(replace(dic["merges"], " " => ""), "(" => ""), ")" => ""), "[" => ""), "]" => ""), ',')) 
    length(bad_merges) == 2n_merges || ERROR(PCDC_CONSTRAINT_ROWMERGE_TAG, "input merges doesn't have matching length")
    merges = Vector{NTuple{2, NAT}}(undef, n_merges)
    for i in eachindex(merges)
        merges[i] = (bad_merges[2i-1], bad_merges[2i])
    end 
    return get_row_merge_constraint(merges)
end

clear_cache(R::RowMergeConstraint) = nothing
known_at_decoding(P::RowMergeConstraint, index::NAT) = true
