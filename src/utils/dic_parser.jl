#=
    Utility function to parse a string into a dictionary.
    The function takes a string s, a prefix string, and a vector of expected field names. 
    It removes the prefix from the string, splits the remaining string into key-value pairs, 
    and constructs a dictionary from these pairs. 
    It also checks that all expected fields are present in the resulting dictionary. 
    
    Warning : the function is currently deprecated and may not handle nested structures correctly.
=#
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
#=
    Utility function to parse a string into a dictionary.
=#
function parse_vector(::Type{T}, s::AbstractString) where {T}
    Vector{T}(parse.(T, split(replace(replace(replace(s, " " => ""), "]" => ""), "[" => ""), ',')))
end