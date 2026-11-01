compress!(u::Vector{BIT}, msb_first::Bool, shift::Bool, rev::Bool, dst::Vector{UInt8}) = begin
    (length(dst) <= length(u) / 8 <= length(dst) + 1) || ERROR("Compress", "$(length(u)) cannot be compressed to length $(length(dst))")
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