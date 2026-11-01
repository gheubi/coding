using CRC
include("compress.jl")
"""
    For small CRC lengths (≤ 15), specs are defined here : 
        https://reveng.sourceforge.io/crc-catalogue/1-15.htm
    For length 16 : http://reveng.sourceforge.net/crc-catalogue/16.htm#crc.cat.crc-16-ccitt-fals
    For length ≥ 17, it is here : http://reveng.sourceforge.net/crc-catalogue/17plus.htm

    For further implementation of missing entries, see e.g. https://users.ece.cmu.edu/~koopman/crc/notes.html
"""
function getSpecs(crclen::Integer)
    if crclen == 3
        return spec(3, 0x7, 0x0, false, false, 0x0, 0x0)#CRC_3_ROHC#spec(3, 0x7, 0x0, false, false, 0x0, 0x0)#CRC_3_ROHC
    elseif crclen == 4
        return spec(4, 0x03, 0x00, false, false, 0x00, 0x0)
    elseif crclen == 5
        return spec(5, 0x09, 0x09, false, false, 0x00, 0x00)
    elseif crclen == 6
        return spec(6, 0x27, 0x3f, false, false, 0x00, 0x0d)
    elseif crclen == 7
        return spec(7, 0x09, 0x00, false, false, 0x00, 0x75)
    elseif crclen == 8
        return spec(8, 0x9b, 0xff, false, false, 0x00, 0xda)
    elseif crclen == 10
        return spec(10, 0x0233, 0x0000, false, false, 0x0000, 0x0199)
    elseif crclen == 11
        return spec(11, 0xe21, 0x000, false, false, 0x0000, 0x0000)
    elseif crclen == 12
        return spec(12, 0x0f13, 0x0fff, false, false, 0x0000, 0x0d4d)
    elseif crclen == 13
        return spec(13, 0x1cf5, 0x0000, false, false, 0x0000, 0x04fa)
    elseif crclen == 14
        return spec(14, 0x0805, 0x0000, false, false, 0x0000, 0x082d)
    elseif crclen == 15
        return spec(15, 0x4599, 0x0000, false, false, 0x0000, 0x059e)
    elseif crclen == 16
        return spec(16, 0xe92f, 0x0000, false, false, 0x0000, 0x0000)
        #return spec(16, 0xf497, 0x0000, false, false, 0x0000, 0x0000)
    elseif crclen == 17
        return spec(17, 0x1685b, 0x00000, false, false, 0x00000, 0x04f03)
    elseif crclen == 21
        return spec(21, 0x102899, 0x000000, false, false, 0x000000, 0x0ed841)
    elseif crclen == 24
        return spec(24, 0x00864cfb, 0x00b704ce, false, false, 0x0000, 0x0021cf02)
    elseif crclen == 30
        return spec(30, 0x2030b9c7, 0x3fffffff, false, false, 0x3fffffff, 0x04c34abf)
    elseif crclen == 31
        return spec(31, 0x04c11db7, 0x7fffffff, false, false, 0x7fffffff, 0x0ce9e46c)
    elseif crclen == 32
        return spec(32, 0x04c11db7, 0xffffffff, false, false, 0xffffffff, 0xcbf43926)
    elseif crclen == 40
        return CRC_40_GSM
    elseif crclen == 64
        return CRC_64
    elseif crclen == 82
        return CRC_82_DARC
    else
        throw(ArgumentError("Unsupported CRC length : $crclen"))
    end
end

getCRC(crclen::NAT) = crc(getSpecs(crclen))

gen_spec(width::NAT, poly::UInt) = spec(width, poly, 0x00, false, false, 0x00, 0x00)
"""
    @brief compute the crc of a given bit vector [b_1, b_2, ..., b_n]

    @param width the width of the remainder of the poly division
    @param poly the polynomial to perform the division. 
        @warning (poly>>width)&1 must be 1 
    @param u the vector to ocmpute the remainder thereofe
"""
function crc_compute!(width::NAT, poly::Vector{BIT}, u::Vector{BIT}, work::Vector{BIT}, dst::Vector{BIT})
    n = length(u)
    n == length(work) || error("work tape is not of the right size (should be $n, is $(length(work)))")
    width == length(dst) || error("destination tape's size does not match width ($(length(dst)) != $width)")
    width + 1 == length(poly) || error("polynomial length doesn't match width (should be $(width+1), is $(length(poly)))")
    poly[1] == 1 || error("msb of poly is not 1")
    width > n && error("length of vector must be at least as big as width ($n < $width)")

    copy!(work, u)
    if width != n
        @inbounds for i in 1:(n-width)
            work[i] == 0 && continue
            for j in 0:width
                work[i+j] ⊻= poly[j+1]
            end
            work[i] == 0 || error("something went wrong during computation (bug)")
        end
    end
    @inbounds for o in 1:width
        dst[o] = work[end-o+1]
    end
    return dst
end

crc_compute(width::NAT, poly::UInt, u::Vector{BIT}, rev::Bool=false, msb_first::Bool=false, shift::Bool=false) =
    crc_compute(width, poly, u, zeros(UInt8, ceil(NAT, length(u)/8)), rev, msb_first, shift)

crc_compute(width::NAT, poly::UInt, u::Vector{BIT}, chunks::Vector{UInt8}, rev::Bool=false, msb_first::Bool=false, shift::Bool=false) =
    crc(gen_spec(width, poly))(compress!(u, msb_first, shift, rev, chunks))