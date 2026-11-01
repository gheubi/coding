using Distributions, LogExpFunctions

struct RandomizedRoundingSC
    llr_mat::Matrix{REAL}
    llr_in::Matrix{Bool}
    depth::NAT
    padded_infobits::Vector{BIT}
end

function get_randomized_rounding_sc(blocklength::NAT)
    ispow2(blocklength) || ERROR(RR_SC_TAG, "$blocklength should be a power of two")
    depth = convert(NAT, log2(blocklength)) + 1
    llr_mat = zeros(REAL, (blocklength, depth))
    llr_in = zeros(Bool, (blocklength, depth))
    padded_infobits = zeros(BIT, blocklength)
    return RandomizedRoundingSC(llr_mat, llr_in, depth, padded_infobits)
end

function randomized_rounding_encode!(
    rrsc::RandomizedRoundingSC,
    init_llrs::Vector{REAL},
    infobits::Vector{BIT},
    infobitsmask::Vector{Bool},
    dynamic_shaping_bit_mask::Vector{Bool},
    dynamic_shaping_function::DynamicConstraint,
    dst::Vector{BIT}
)
    n = length(init_llrs)

    #checl args
    @check ispow2(n) || ERROR(RR_SC_TAG, "$N is not a power of 2")
    @check length(infobits) == count(infobitsmask) || ERROR(RR_SC_TAG, "length don't match")

    #reinitializing datastruct
    llr_mat = rrsc.llr_mat
    llr_in = rrsc.llr_in
    depth = rrsc.depth
    padded_infobits = rrsc.padded_infobits
    fill!(llr_mat, 0)
    fill!(llr_in, false)
    fill!(padded_infobits, BIT(0))


    #at depth of the tree are the raw LLR of the channel for each output y[i]
    cnt = 1
    @inbounds for i in 1:n
        llr_mat[i, depth] = init_llrs[i]
        llr_in[i, depth] = true
        #padded input vector
        if infobitsmask[i]
            padded_infobits[i] = infobits[cnt]
            cnt += 1
        end
    end

    #small recursive function to get the LLR at a given node (i,d) as a function of U
    #we're using memoization here, i.e when an entry of the llr_mat is computed once it will be only once
    function get(u::Vector{BIT}, i::NAT, d::NAT)::REAL
        shift = n  ÷ (1 << d)
        j = length(u)
        if !llr_in[i, d]
            l2 = get(u[1:2:j-1] .⊻ u[2:2:j], i, d+1)
            newi = i + ((i-1) % 2shift < shift ? shift : -shift)
            l1 = get(u[2:2:j], newi , d+1)
            if iseven(j)
                llr_mat[i, d] = logaddexp(l1+l2,0) - logaddexp(l1, l2)
            else
                llr_mat[i, d] = l2 + (u[j] == 0 ? l1 : -l1)
            end
            llr_in[i, d] = true
        end
        llr_mat[i, d]
    end
    #SC decoder
    for i in 1:n
        #the tree works in bitreversed order
        l = get(dst[1:i-1], bit_rev_op(i, depth-1), 1)
        if infobitsmask[i]
            dst[i] = padded_infobits[i]
        elseif  dynamic_shaping_bit_mask[i]
            dst[i] = dynamic_shaping_function(dst, infobitsmask, i)
        else
            p1 = 1/(1 + exp(l))
            dst[i] = BIT(rand(Bernoulli(p1)))
        end
    end
    @check all(llr_in) || ERROR(RR_SC_TAG, "Somethin bad happend")
    return dst
end