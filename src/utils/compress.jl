#=
    Compress a vector of bits u into a vector of bytes dst.
    The function takes the following parameters:
    - u: a vector of bits (0 or 1) to be compressed
    - msb_first: a boolean indicating whether to store the most significant bit first (true) or least significant bit first (false)
    - shift: a boolean indicating whether to shift the bits or not /
        (shifting means that the last bits will be shifted at the msb of the last byte, non-shifting means that the last bits will placed the lsb of the last byte)
    - rev: a boolean indicating whether to reverse the order of the bits or not 
    - dst: a preallocated vector of bytes to store the compressed output
    The function returns the compressed vector of bytes dst.
=#
compress!(u::Vector{BIT}, msb_first::Bool, shift::Bool, rev::Bool, dst::Vector{UInt8}) = begin
    all(0 .<= u .<= 1) || ERROR(UTILS_COMPRESS_TAG, "$u is not a valid bit vector")
    (length(dst) <= length(u) / 8 <= length(dst) + 1) || ERROR(UTILS_COMPRESS_TAG, "$(length(u)) cannot be compressed to length $(length(dst))")
    chunks = dst
    fill!(chunks, 0)
    k = length(u)
    first_unread_index = rev ? k : 1
    @inbounds for i in eachindex(chunks)
        chunk = chunks[i]
        range = !rev ? 
            msb_first ? (first_unread_index:(first_unread_index+7)) : ((first_unread_index+7):-1:first_unread_index) :
            msb_first ? ((first_unread_index-7):first_unread_index) : (first_unread_index:-1:(first_unread_index-7))
        for j in range 
            !shift && (j > k) && continue
            chunk <<= 1
            shift && (j > k) && continue
            chunk |= 1&u[j]
        end
        first_unread_index += rev ? -8 : 8
        chunks[i] = chunk
    end
    return chunks
end

compress(u::Vector{BIT}, msb_first::Bool, shift::Bool, rev::Bool) = compress!(u, msb_first, shift, rev, zeros(BIT, ceil(NAT, length(u)/8)))