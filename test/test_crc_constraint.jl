include("../conf.jl")

function crc_compute!(width::NAT, poly::Vector{BIT}, u::Vector{BIT}, work::Vector{BIT}, dst::Vector{BIT})
    copy!(work, u)
    @inbounds for i in 1:(length(u)-width)
        work[i] == 0 && continue
        for j in 0:width
            work[i+j] ⊻= poly[j+1]
        end
        work[i] == 0 || ERROR(crc_c_tag, "something went wrong during crc-computation (bug)")
    end
    @inbounds for o in 1:width
        dst[o] = work[end-o+1]
    end
    return dst
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
        ERROR(crc_c_tag, "default polynomial not added for width $width")
    end
end



for width in [3, 4, 5, 11, 16]
    dict = Dict{UInt, Vector{BIT}}()
    poly = get_default_poly(width)
    work = zeros(BIT, 32)
    dst = zeros(BIT, width)
    for _ in 1:1000
        vec = BIT.(rand(Bool, 32))
        h = hash(vec)
        crc = get!(dict, h) do
            crc_compute!(width, poly, vec, work, dst)
        end
        crc = copy(crc)
        vec2 = BIT.(rand(Bool, 32))
        h2 = hash(vec2)
        crc2 = get!(dict, h2) do
            crc_compute!(width, poly, vec2, work, dst)
        end
        crc3 = get!(dict, h) do
            crc_compute!(width, poly, vec, work, dst)
        end
        if crc != crc3
            print("$crc != $crc3")
        end
    end
end