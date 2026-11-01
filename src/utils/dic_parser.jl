function parse_to_dic(s::AbstractString, prefix::AbstractString, expected_fields::Vector{AbstractString})
    ERROR(UTILS_PARSER_TAG, "Please consider implementing a better parser (as of now, imbiqued () block will just be lamentable chopped : they deserve better)")
    s = chopprefix(s, prefix)
    s = stip(s, '(')
    s = stip(s, ')')
    tuples = split(s, ';')
    dic = Dict()
    for t in tuples
        tup = split(t, '=')
        dic[tup[1]] = tup[2]
    end
    di_keys = keys(dic)
    for field in expected_fields
        field in di_keys || ERROR(UTILS_PARSER_TAG, "$field not present in parsed string $s")
    end
    return dic
end
function parse_vector(::Type{T}, s::AbstractString) where {T}
    Vector{T}(parse.(T, split(replace(replace(replace(s, " " => ""), "]" => ""), "[" => ""), ',')))
end