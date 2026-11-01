include("fetch_rowmerge.jl")

using RowEchelon

# determines wether i ⪯ j according to the partial order on bit-channels
⪯(i::NAT, j::NAT) = begin
    function binarydom(i0, j0)
        m = max(i0, j0)
        while m > 0
            il = i0&1
            jl = j0&1
            i0 >>>= 1
            j0 >>>= 1
            m >>>= 1
            il <= jl || return false
        end
        return true
    end
    function leftswap(i0, j0)
        merged = i0⊻j0
        weight(merged) == 2 || return false
        while 1&merged == 0
            merged >>>= 1
            i0 >>>= 1
            j0 >>>= 1
        end
        i0 == 1 && j0 == 2 && return true
        return false
    end
    binarydom(i-1, j-1) || leftswap(i-1, j-1)
end

# a rate profile q is said to be decreasing if
# for all i ∈ q and j ∈ [1, n] : i ≤ j ⟹ j ∉ q
is_decreasing(q::Vector{NAT}, N::NAT) = begin
    for i ∈ q,j ∈ 1:N
        if i ⪯ j && j ∉ q
            return false
        end      
    end
    return true
end

# a rate profile q is said to be decreasing if
# for all i ∈ q and j ∈ [1, n] : i ≤ j ⟹ j ∉ q
is_decreasing(paddedp::Vector{Bool}) = begin
    N = length(paddedp)
    @inbounds for i in eachindex(paddedp)
        if paddedp[i]
            for j in 1:N
                if i ⪯ j
                    if !paddedp[i]
                        return false
                    end
                end
            end
        end
    end
    return true
end

# returns the weigth of the row of the polar trasform at index i
# @warning this function use 0-based indexing
function polar_row_weight_0(i0::NAT)
    1 << weight(i0)
end

# returns the weigth of the merges rows at index i and j
# @waring this function use 0-based indexing
function merged_polar_row_weight_0(i0::NAT, j0::NAT)
    #=
        See theorem 1 of Gelincik : "Achieving PAC Code Performance 
        with SCL decoding without extra computational complexity"
    =#
    polar_row_weight_0(i0) + polar_row_weight_0(j0) - 2*polar_row_weight_0(i0&j0)
end

# finds the frozen rows that increasing the weigths
function increasing_weight_frozen_rows!(padded_info_set::Vector{Bool}, i::NAT, dst1::Vector{Bool}, dst2::Vector{Bool})
    N = length(padded_info_set)
    N0 = N-1
    i0 = i-1
    base_w = polar_row_weight_0(i0)
    @inbounds for f0 in 0:N0
        if f0 > i0 && !padded_info_set[f0+1]
            b = merged_polar_row_weight_0(i0, f0) > base_w
            dst1[f0+1] = b
            dst2[f0+1] = !b
        else 
            dst1[f0+1] = false
            dst2[f0+1] = false
        end
    end
    return dst1, dst2
end
increasing_weight_frozen_rows(padded_info_set::Vector{Bool}, i::NAT) = 
    increasing_weight_frozen_rows!(padded_info_set, i, zeros(Bool, length(padded_info_set)), zeros(Bool, length(padded_info_set)))
     

function min_dist(padded_info_set::Vector{Bool}, merging)
    if is_decreasing(padded_info_set)
        return min_dist(padded_info_set)
    else 
        error("Implement me!")
    end
end

function min_dist(padded_info_set::Vector{Bool})
    acc = 0
    inited = false
    @inbounds for i ∈ 1:length(padded_info_set)
        if padded_info_set[i]
            w = polar_row_weight_0(i-1)
            if !inited
                acc = w
                inited = true
            end
            if w < acc
                acc = w 
            end
        end
    end
    return acc
end
function padded_weight_index!(N::NAT, w::NAT, dst::Vector{Bool})
    @inbounds for i ∈ 1:N
        dst[i] = polar_row_weight_0(i-1) == w
    end
    return dst
end
padded_weight_index(N::NAT, w::NAT) = padded_weight_index!(N, w, zeros(Bool, N))
not_zero(i0, N0) = (N0)⊻(i0)
in_ki0(j0, i0, N0) = begin
    #N >= j > i && mergedPolarRowWeigth(i,j) == polarRowWeigth(i)
    weight(not_zero(i0, N0)&j0) == 1
end
not_in_ki0(j0,i0,N0) = begin
    #N >= j > i && !(mergedPolarRowWeigth(i,j) == polarRowWeigth(i))
    weight(not_zero(i0, N0)&(j0)) > 1
end
function max_padded(padded_vector::Vector{Bool})
    i = length(padded_vector)
    while i > 0 && !padded_vector[i]
        i -= 1
    end
    (i == 0) && error("No max over an empty set")
    return i
end
function min_weight_number(padded_info_set::Vector{Bool}, min_weight::NAT)
    #=
        This assumes the infoset I satisfies the partial order property.
        We here rely on the paper by Rowshan : On the formation of min-weight
        Codewords of Polar/PAC codes and its applications.
    =#
    N = length(padded_info_set)
    N0 = N-1
    n = convert(NAT, log2(N))
    n0 = n-1
    min_weight_n = 0
    @inbounds for i in eachindex(padded_info_set)
        i0 = i-1
        if padded_info_set[i] && polar_row_weight_0(i0) == min_weight
            #println("Found $i ∈ I_wmin : $(bitstring(i-1))")
            ki_length = weight(not_zero(i0, N0))
            for k0 in 0:n0
                if (i0>>k0)&1 == 1
                    #println("Found that $k ∈ supp($(i-1))")
                    for l in (k0+1):n0
                        #println("for $l : ¬i_l = $(((i0>>l)&1)⊻1)")
                        ki_length += ((i0>>l)&1)⊻1
                    end
                end 
            end
            #println("|K_i| = $ki_length")
            min_weight_n += (1 << ki_length)
            #println("New min weightN = $min_weight_n")
        end
    end
    return min_weight_n
end
function min_weight_number_3(padded_info_set::Vector{Bool}, min_weight::NAT)
    length_n = 0
    N = length(padded_info_set)
    N0 = N-1
    @inbounds for i ∈ eachindex(padded_info_set)
        i0 = i-1
        if padded_info_set[i] && polar_row_weight_0(i0) == min_weight
            length_k = 0
            for j0 ∈ i0+1:N0
                if polar_row_weight_0(j0) >= merged_polar_row_weight_0(i0,j0) == min_weight
                    length_k += 1
                end
            end
            length_n += (1 << length_k)
        end
    end
    return length_n
end
function min_weight_number(padded_info_set::Vector{Bool}, merging::Vector{NTuple{2, NAT}}, wmin::NAT)
    N = length(padded_info_set)
    N0 = N-1
    μ0(i0, j0, k0) = begin
        return (not_zero(i0, N0)&(j0|k0))|(j0&k0)
    end
    function update_message(i0,j0,u)
        for k0 ∈ (i0+1):(j0-1) 
            if u[k0+1] == 1 && (not_zero(i0,N0)&(j0)&(k0)) == 0
                u[μ0(i0,j0,k0)+1] ⊻= 1 
            end
        end
        u[j0+1] = 1
        return u
    end
    xor(i0,k0, u) = begin
        xacc = 0
        for h0 ∈ i0:(k0-1)
            xacc += u[h0+1]*T[h0+1,k0+1]
        end
        return xacc % 2
    end
    fistar0(i0) = begin
        j0 = N-1
        while j0 > i0
            if not_in_ki0(j0,i0,N0) && !padded_info_set[j0+1]
                return j0
            end
            j0 -= 1
        end
        return i0
    end
    function enumarte_subtree(i0, j0, u)
        #println("Here we go :D ($i, $j, $u)")
        awmin = 0
        u = update_message(i0,j0,copy(u))
        for k0 ∈ (j0+1):fistar0(i0)
            if padded_info_set[k0+1] && in_ki0(k0,i0,N0)
                awmin += enumarte_subtree(i0, k0, u)
            elseif !padded_info_set[k0+1] && u[k0+1] != xor(i0,k0,u)
                if in_ki0(k0,i0,N0)
                    u = update_message(i0,k0,u)
                else
                    return awmin
                end
            end 
        end
        return awmin + 1
    end
    #print("Compute matrix representation...")
    T = zeros(BIT, (N,N))
    @inbounds for i ∈ eachindex(padded_info_set)
        T[i, i] = padded_info_set[i]
    end
    for t ∈ merging
        T[t[1], t[2]] = 1
    end 
    T[padded_info_set, :] .= rref!(T[padded_info_set, :])
    #println("done")
    nweight = 0
    for i ∈ eachindex(padded_info_set)
        i0 = i-1
        if padded_info_set[i] && polar_row_weight_0(i0) == wmin
            #println("Found $i ∈ I_wmin")
            f0 = fistar0(i0)
            #println("Compute fstar : $f")
            expo = count(j0 -> begin
                padded_info_set[j0+1] && in_ki0(j0, i0, N0)
            end, f0:N0)
            #println("Computed exp $expo")
            c = enumarte_subtree(i0, i0, zeros(BIT, N))
            nweight += c*(1<<expo)
        end
    end
    return nweight
end

#=
    Compute optimal merges for the info set padded_info_set by finding pairs of indices (r, m) such that the weigth 
    spectrum of the resulting code is improved. The algorithm is forbidden to select merges from the forbidden set

    @note this function should never be directly called, please prefer its alternative signature who saves compuation results
    @note see Zunker et al, "Row-Merged Polar Codes: Analysis, Design, and Decoder Implementation" for the origin of the algorithm
    
    @warning this function doesn't give any optimality guarrantee when the forbiden_set is not empty, or if the padded_info_set
        is not decreasing (in the sense of is_decreasing())

    @return a set of merges R = {(r, m) | r in padded_info_set}
=#
function merge_rows(padded_info_set::Vector{Bool}, forbiden_set::Vector{NAT})
   N = length(padded_info_set)

    D = zeros(Bool, N)
    opt_d = zeros(Bool, N)

    R = Vector{NTuple{2,NAT}}([])
    opt_r = R
    
    i_trans = zeros(Bool, N)
    f_star = zeros(Bool, N)
    f_loser = zeros(Bool, N)
    f_trans = f_star


    dmin = min_dist(padded_info_set, R)
    dmin_trans= 0
    opt_dmin = 0

    admin = min_weight_number(padded_info_set, R, dmin)
    admin_trans = 0
    opt_admin = 0

    increase_weight = true
    
    while true
        opt_dmin = dmin
        opt_admin = admin
        @LOG ROWMERGE_TAG "On main loop with increase_weight=$increase_weight : opt_dmin=$opt_dmin, opdAdmin = $opt_admin"
        padded_weight_index!(N, opt_dmin, i_trans)
        for i ∈ eachindex(i_trans)
            if i_trans[i]
                increasing_weight_frozen_rows!(padded_info_set, i, f_star, f_loser)
                all(j -> !j, f_star) && continue

                if increase_weight
                    f_star .&= .!D
                    f_trans = f_star
                else
                    f_loser .|= f_star
                    f_loser .&= .!D
                    f_trans = f_loser
                end
                c = 0
                for b in eachindex(f_trans)
                    if c <= 2 && f_trans[b] && b ∉ forbiden_set
                        b <= i && ERROR(ROWMERGE_TAG, "Should consider higher indices for merging rows")
                        c += 1
                        r_trans = copy(R)
                        t::NTuple{2, NAT} = NTuple{2, NAT}((i,b))
                        push!(r_trans, t)
                        @LOGD row_merge_tag "   Found merge candidate $t"
                        dmin_trans = min_dist(padded_info_set, r_trans)
                        if dmin_trans >= opt_dmin
                            @LOGD ROWMERGE_TAG "     Computing admin for a merge at index $b..."
                            admin_trans = min_weight_number(padded_info_set, r_trans, dmin_trans)
                            @LOG ROWMERGE_TAG "w=$dmin_trans admin=$admin_trans"
                            if dmin_trans > opt_dmin || admin_trans < opt_admin
                                @LOG ROWMERGE_TAG "Merges improve results : new opt_dmin=$dmin_trans, new opt_admin=$admin_trans"
                                opt_dmin = dmin_trans
                                opt_admin = admin_trans
                                copy!(opt_d, D)
                                opt_d[b] = true
                                opt_r = r_trans
                            end
                        end
                    end
                end
            end
        end
        if opt_dmin > dmin || opt_admin < admin
            copy!(D,opt_d)
            R = opt_r
            dmin = opt_dmin
            admin = opt_admin
        elseif increase_weight
            increase_weight = false
        else
            @LOG ROWMERGE_TAG "Found $(length(R)) merges"
            return R
        end
    end
end
#=
    Compute optimal merges for the info set padded_info_set by finding pairs of indices (r, m) such that the weigth 
    spectrum of the resulting code is improved. The algorithm is forbidden to select merges from the forbidden set

    @note see Zunker et al, "Row-Merged Polar Codes: Analysis, Design, and Decoder Implementation" for the origin of the algorithm
    
    @warning this function doesn't give any optimality guarrantee when the forbiden_set is not empty, or if the padded_info_set
        is not decreasing (in the sense of is_decreasing())

    @return a set of merges R = {(r, m) | r in padded_info_set}
=#
function merge_rows(info_set::Vector{NAT}, N::NAT, forbiden_set::Vector{NAT})
    id =  merges_already_present(N, info_set, forbiden_set)
    if id > 0 
        @LOG ROWMERGE_TAG "merges already computed"
        return merges_fetch(id)
    end
    padded_info_set = zeros(Bool, N)
    @inbounds for i in info_set
        padded_info_set[i] = true
    end
    merges = merge_rows(padded_info_set, forbiden_set)
    dump_merges(N, info_set, forbiden_set, merges)
    return merges
end
#=
    Compute optimal merges for the info set padded_info_set by finding pairs of indices (r, m) such that the weigth 
    spectrum of the resulting code is improved.

    @note see Zunker et al, "Row-Merged Polar Codes: Analysis, Design, and Decoder Implementation" for the origin of the algorithm
    
    @warning this function doesn't give any optimality guarrantee if the padded_info_set
        is not decreasing (in the sense of is_decreasing())

    @return a set of merges R = {(r, m) | r in padded_info_set}
=#
merge_rows(info_set::Vector{NAT}, N::NAT) = merge_rows(info_set, N, zeros(NAT, 0))

# Find the message length K of a Reed-Muller code RM(r,n)
function rm_k(r::NAT, n::NAT)
    acc = 0
    for l ∈ 0:r
        acc += binomial(n,l)
    end
    return acc
end
# Find the smallest r of a Reed-Muller code RM(r, n) with a messagelength >= K
function find_least_r(n::NAT, K::NAT)
    #generated by chatGPT
    S = 0
    c = 1
    for r in 0:n
        S += c
        S > K && return r
        r < n && (c = div(c*(n-r), r+1))
    end
    return -1
end

# Select the K indices of an information set that have the highest reliabilites and such that the
# weight of the rows of the polar transform is optimized. 
function target_next_rm_code(N::NAT, K::NAT, reliabilities::Vector{REAL}, forbiden_set::Vector{NAT})
    ranking = sort(reverse(1:N), by=i->reliabilities[i])
    n = convert(NAT, log2(N))
    k = K
    while k <= N
        r = find_least_r(n, k)
        r < 0 && break 
        threshold = (1 << (n-r))
        rm_info_set = filter(i -> polar_row_weight_0(i-1) >= threshold && i ∉ forbiden_set, ranking)
        length(rm_info_set) >= K && return NAT.(rm_info_set[1:K])
        k += 1
    end
    ERROR(ROWMERGE_TAG, "Couldn't find any merges for N=$N, K=$K")
end
target_next_rm_code(N::NAT, K::NAT, reliabilites::Vector{REAL}) = target_next_rm_code(N, K, reliabilites, zeros(NAT, 0))
