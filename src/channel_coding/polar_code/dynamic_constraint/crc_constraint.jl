struct CrcConstraint <: DynamicConstraint
    _width::NAT
    _poly::Vector{BIT}

    _work_tapes::Vector{Vector{BIT}}

    _index_perm::Vector{NAT}
    _visible_at_decoding_indices::Vector{Bool}

    _cache::Vector{Vector{Tuple{UInt, Vector{BIT}}}}
    _version_number::NAT
end

function get_default_poly(width::NAT)
    if width == 3 
        return BIT.([1, 0, 1, 1])
    elseif width == 4
        return BIT.([1, 0, 0, 1, 1])
    elseif width == 5
        return BIT.([1, 0, 1, 0, 0, 1])
    elseif width == 11
        return BIT.([1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    elseif width == 16
        return BIT.([1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1])
    else
        ERROR(PCDC_CONSTRAINT_CRC_TAG, "default polynomial not added for width $width")
    end
end

function get_crc_constraint(width::NAT, poly::Vector{BIT}, K::NAT, indexes::Vector{NAT}, visibles::Vector{Bool},version_number::NAT=1)
    @LOGV PCDC_CONSTRAINT_CRC_TAG "generating constraints"
    length(indexes) == width || ERROR(PCDC_CONSTRAINT_CRC_TAG, "indexes don't have matching size")
    length(unique(indexes)) == length(indexes) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "indexes contain duplicates")
    length(visibles) == width || ERROR(PCDC_CONSTRAINT_CRC_TAG, "visibles indices don't have matching size")
    width+1 == length(poly) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "polynomial length doesn't match width (should be $(width+1), is $(length(poly)))")
    poly[1] == 1 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "polynomial not valid : msb (index 1) must 1")
    1 <= version_number <= 2 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "unknown version number $(version_number)")

    work_tapes = Vector{Vector{BIT}}(undef, NDEFAULT)
    dic = Vector{Vector{Tuple{UInt, Vector{BIT}}}}(undef, NDEFAULT)
    @inbounds for i in eachindex(work_tapes)
        if version_number == 1
            work_tapes[i] = zeros(BIT, K)
        elseif version_number == 2
            work_tapes[i] = zeros(BIT, K+width)
        end
        dic[i] = Vector{Tuple{UInt, Vector{BIT}}}(undef, 32)
        @inbounds for j in 1:32
            dic[i][j] = (UInt(j), zeros(BIT, width))
        end
    end 
    return CrcConstraint(width, poly, work_tapes, indexes, visibles, dic, version_number)
end

get_crc_constraint(width::NAT, poly::Vector{BIT}, K::NAT, indexes::Vector{NAT},version_number::NAT=1) = get_crc_constraint(width, poly, K, indexes, zeros(Bool, width),version_number)

get_crc_constraint(width::NAT, K::NAT, indexes::Vector{NAT},version_number::NAT=1) = get_crc_constraint(width, get_default_poly(width), K, indexes,version_number)

function compress!(u::Vector{BIT}, select_index::Vector{Bool}, dst::Vector{BIT})
    (count(select_index) == length(dst) && length(u) == length(select_index)) || 
        ERROR(PCDC_CONSTRAINT_CRC_TAG, "cannot compress vector ($(count(select_index)) selected index for a dest of length $(length(dst)))")

    dst_counter = 1
    @inbounds for i in eachindex(u)
        if select_index[i]
            dst[dst_counter] = u[i]
            dst_counter += 1
        end
    end
    return dst
end

"""
    @brief compute the crc of a given bit vector [b_1, b_2, ..., b_n]

    @param width the width of the remainder of the poly division
    @param poly the polynomial to perform the division. 
    @param u the vector to ocmpute the remainder thereofe
"""
function crc_compute_v1!(width::NAT, poly::Vector{BIT}, work::Vector{BIT}, dst::Vector{BIT})
    n = length(work)
    @check width == length(dst) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "destination tape's size does not match width ($(length(dst)) != $width)")
    @check width+1 == length(poly) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "polynomial length doesn't match width (should be $(width+1), is $(length(poly)))")
    @check poly[1] == 1 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "msb of poly is not 1")
    @check width > n && ERROR(PCDC_CONSTRAINT_CRC_TAG, "length of vector must be at least as big as width ($n < $width)")

    if width != n
        @LOGV PCDC_CONSTRAINT_CRC_TAG "xoring a bunch of times"
        @inbounds for i in 1:(n-width)
            work[i] == 0 && continue
            for j in 0:width
                work[i+j] ⊻= poly[j+1]
            end
            work[i] == 0 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "something went wrong during crc-computation (bug, please solve me)")
        end
    end
    @inbounds for o in 1:width
        dst[o] = work[end-o+1]
    end
    @LOGV PCDC_CONSTRAINT_CRC_TAG "computed crc-digest : $(print_vector(dst))"
    return dst
end

function crc_compute!(width::NAT, poly::Vector{BIT}, work::Vector{BIT}, dst::Vector{BIT})
    extended = length(work)
    n = extended - width
    @check width == length(dst) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "destination tape's size does not match width ($(length(dst)) != $width)")
    @check width+1 == length(poly) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "polynomial length doesn't match width (should be $(width+1), is $(length(poly)))")
    @check poly[1] == 1 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "msb of poly is not 1")
    @check all(i -> work[i] == 0, n+1:extended) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "m[x] was not properly shifted")

    @LOGV PCDC_CONSTRAINT_CRC_TAG "xoring a bunch of times"
    @inbounds for i in 1:n
        work[i] == 0 && continue
        for j in 0:width
            work[i+j] ⊻= poly[j+1]
        end
        work[i] == 0 || ERROR(PCDC_CONSTRAINT_CRC_TAG, "something went wrong during crc-computation (bug, please solve me)")
    end
    @inbounds for o in 1:width
        dst[o] = work[n+o]
    end
    @LOGV PCDC_CONSTRAINT_CRC_TAG "computed crc-digest : $(print_vector(dst))"
    return dst
end

function translate_index(C::CrcConstraint, index::NAT)
    @inbounds for i in eachindex(C._index_perm)
        if C._index_perm[i] == index
            return i
        end
    end
    ERROR(PCDC_CONSTRAINT_CRC_TAG, "$index is not a CRC-index")
end

function get_or_compute(C::CrcConstraint, h::UInt, work_tape::Vector{BIT})
    tid = default_tid()
    v = C._cache[tid]
    index = (h % 32) + 1
    h2, v2 = v[index]
    h2 == h && return v2
    res =  C._version_number == 1 ? 
        crc_compute_v1!(C._width, C._poly, work_tape, v2) : 
        crc_compute!(C._width, C._poly, work_tape, v2)
    v[index] = (h, res)
    return res 
end


(C::CrcConstraint)(u::Vector{BIT}, select_index::Vector{Bool}, index::NAT) = begin
    @LOGV PCDC_CONSTRAINT_CRC_TAG "Computing dynamic bit at $index"
    @check count(select_index) == length(C._work_tapes[1]) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "length of selected indices doesn't match work tapes")
    @check length(u) == length(select_index) || ERROR(PCDC_CONSTRAINT_CRC_TAG, "length of input vector u doesn't match that of selected index")

    index = translate_index(C, index)
    tid = default_tid()
    @LOGV PCDC_CONSTRAINT_CRC_TAG "compressing input"
    work_tape = C._work_tapes[tid]
    C._version_number >= 2 && fill!(work_tape, BIT(0))
    work_tape = compress!(u, select_index, work_tape)
    h = hash(work_tape)
    @LOGV PCDC_CONSTRAINT_CRC_TAG "computing $index-th bit of crc with input hash $h, looking up in dictionnary"
    res = get_or_compute(C, h, work_tape)
    return res[index]
end

known_at_decoding(C::CrcConstraint, index::NAT) = begin
    index = translate_index(C, index)
    return C._visible_at_decoding_indices[index]
end

clear_cache(C::CrcConstraint) = begin
    v = C._cache[default_tid()]
    @inbounds for i in eachindex(v)
        h, v2 = v[i]
        v[i] = (UInt(i), v2)
    end
end


Base.show(io::IO, C::CrcConstraint) = begin
    C._version_number >= 2 ?  print(io, "CRC(width=$(C._width);poly=$(Int.(C._poly));input_length=$(length(C._work_tapes[1]));index_perm=$(C._index_perm);known_at_decoding=$(NAT.(C._visible_at_decoding_indices));v=$(C._version_number))") :
        print(io, "CRC(width=$(C._width);poly=$(Int.(C._poly));input_length=$(length(C._work_tapes[1]));index_perm=$(C._index_perm);known_at_decoding=$(NAT.(C._visible_at_decoding_indices)))")
end
Base.parse(::Type{T}, s::AbstractString) where {T <: CrcConstraint} = begin
    startswith(s, "CRC(") || ERROR(PCDC_CONSTRAINT_CRC_TAG, "input string cannot be parsed")

    dic = parse_to_dic(s, "CRC", ["width", "poly", "input_length", "index_perm", "visible_at_decoding_indices"])
    width = parse(NAT, dic["width"])
    poly = parse_vector(BIT, dic["poly"])
    input_length = parse(Int, dic["input_length"])
    index_perm = parse_vector(NAT, dic["index_perm"])
    visible_at_decoding_indices = parse_vector(NAT, dic["visible_at_decoding_indices"])
    v = "v" in keys(dic) ? parse(NAT,dic["v"]) : 1

    return get_crc_constraint(width, poly, input_length, index_perm, visible_at_decoding_indices,v)
end

