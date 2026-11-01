include("fetch_rowmerge.jl")

using RowEchelon


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

isDecreasing(q::Vector{NAT}, N::NAT) = begin
    for i ∈ q,j ∈ 1:N
        if i ⪯ j && j ∉ q
            return false
        end      
    end
    return true
end
isDecreasing(paddedp::Vector{Bool}) = begin
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

function polarRowWeigth0(i0::NAT)
    1 << weight(i0)
end
function mergedPolarRowWeigth0(i0::NAT, j0::NAT)
    """
        See theorem 1 of Gelincik : "Achieving PAC Code Performance 
        with SCL decoding without extra computational complexity"
    """
    polarRowWeigth0(i0) + polarRowWeigth0(j0) - 2*polarRowWeigth0(i0&j0)
end

function increasingWeightFrozenRows!(paddedInfoSet::Vector{Bool}, i::NAT, dst1::Vector{Bool}, dst2::Vector{Bool})
    N = length(paddedInfoSet)
    N0 = N-1
    i0 = i-1
    baseW = polarRowWeigth0(i0)
    @inbounds for f0 in 0:N0
        if f0 > i0 && !paddedInfoSet[f0+1]
            b = mergedPolarRowWeigth0(i0, f0) > baseW
            dst1[f0+1] = b
            dst2[f0+1] = !b
        else 
            dst1[f0+1] = false
            dst2[f0+1] = false
        end
    end
    return dst1, dst2
end
increasingWeightFrozenRows(paddedInfoSet::Vector{Bool}, i::NAT) = 
    increasingWeightFrozenRows!(paddedInfoSet, i, zeros(Bool, length(paddedInfoSet)), zeros(Bool, length(paddedInfoSet)))
     

function minDist(paddedInfoSet::Vector{Bool}, merging)
    if isDecreasing(paddedInfoSet)
        return minDist(paddedInfoSet)
    else 
        error("Implement me!")
    end
end

function minDist(paddedInfoSet::Vector{Bool})
    acc = 0
    inited = false
    @inbounds for i ∈ 1:length(paddedInfoSet)
        if paddedInfoSet[i]
            w = polarRowWeigth0(i-1)
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
function paddedWeightIndex!(N::NAT, w::NAT, dst::Vector{Bool})
    @inbounds for i ∈ 1:N
        dst[i] = polarRowWeigth0(i-1) == w
    end
    return dst
end
paddedWeightIndex(N::NAT, w::NAT) = paddedWeightIndex!(N, w, zeros(Bool, N))
not0(i0, N0) = (N0)⊻(i0)
inKi0(j0, i0, N0) = begin
    #N >= j > i && mergedPolarRowWeigth(i,j) == polarRowWeigth(i)
    weight(not0(i0, N0)&j0) == 1
end
notInKi0(j0,i0,N0) = begin
    #N >= j > i && !(mergedPolarRowWeigth(i,j) == polarRowWeigth(i))
    weight(not0(i0, N0)&(j0)) > 1
end
function maxPadded(paddedVector::Vector{Bool})
    i = length(paddedVector)
    while i > 0 && !paddedVector[i]
        i -= 1
    end
    (i == 0) && error("No max over an empty set")
    return i
end
function minWeightNumber(paddedInfoSet::Vector{Bool}, minWeight::NAT)
    """
        This assumes the infoset I satisfies the partial order property.
        We here rely on the paper by Rowshan : On the formation of min-weight
        Codewords of Polar/PAC codes and its applications.
    """
    N = length(paddedInfoSet)
    N0 = N-1
    n = convert(NAT, log2(N))
    n0 = n-1
    minWeightN = 0
    @inbounds for i in eachindex(paddedInfoSet)
        i0 = i-1
        if paddedInfoSet[i] && polarRowWeigth0(i0) == minWeight
            #println("Found $i ∈ I_wmin : $(bitstring(i-1))")
            KiLength = weight(not0(i0, N0))
            for k0 in 0:n0
                if (i0>>k0)&1 == 1
                    #println("Found that $k ∈ supp($(i-1))")
                    for l in (k0+1):n0
                        #println("for $l : ¬i_l = $(((i0>>l)&1)⊻1)")
                        KiLength += ((i0>>l)&1)⊻1
                    end
                end 
            end
            #println("|K_i| = $KiLength")
            minWeightN += (1 << KiLength)
            #println("New min weightN = $minWeightN")
        end
    end
    return minWeightN
end
function minWeightNumber3(paddedInfoSet::Vector{Bool}, minWeight::NAT)
    lengthN = 0
    N = length(paddedInfoSet)
    N0 = N-1
    @inbounds for i ∈ eachindex(paddedInfoSet)
        i0 = i-1
        if paddedInfoSet[i] && polarRowWeigth0(i0) == minWeight
            lenghtK = 0
            for j0 ∈ i0+1:N0
                if polarRowWeigth0(j0) >= mergedPolarRowWeigth0(i0,j0) == minWeight
                    lenghtK += 1
                end
            end
            lengthN += (1 << lenghtK)
        end
    end
    return lengthN
end
function minWeightNumber(paddedInfoSet::Vector{Bool}, merging::Vector{NTuple{2, NAT}}, wmin::NAT)
    N = length(paddedInfoSet)
    N0 = N-1
    μ0(i0, j0, k0) = begin
        return (not0(i0, N0)&(j0|k0))|(j0&k0)
    end
    function updateMessage(i0,j0,u)
        for k0 ∈ (i0+1):(j0-1) 
            if u[k0+1] == 1 && (not0(i0,N0)&(j0)&(k0)) == 0
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
            if notInKi0(j0,i0,N0) && !paddedInfoSet[j0+1]
                return j0
            end
            j0 -= 1
        end
        return i0
    end
    function enumarteSubtree(i0, j0, u)
        #println("Here we go :D ($i, $j, $u)")
        awmin = 0
        u = updateMessage(i0,j0,copy(u))
        for k0 ∈ (j0+1):fistar0(i0)
            if paddedInfoSet[k0+1] && inKi0(k0,i0,N0)
                awmin += enumarteSubtree(i0, k0, u)
            elseif !paddedInfoSet[k0+1] && u[k0+1] != xor(i0,k0,u)
                if inKi0(k0,i0,N0)
                    u = updateMessage(i0,k0,u)
                else
                    return awmin
                end
            end 
        end
        return awmin + 1
    end
    #print("Compute matrix representation...")
    T = zeros(BIT, (N,N))
    @inbounds for i ∈ eachindex(paddedInfoSet)
        T[i, i] = paddedInfoSet[i]
    end
    for t ∈ merging
        T[t[1], t[2]] = 1
    end 
    T[paddedInfoSet, :] .= rref!(T[paddedInfoSet, :])
    #println("done")
    nweight = 0
    for i ∈ eachindex(paddedInfoSet)
        i0 = i-1
        if paddedInfoSet[i] && polarRowWeigth0(i0) == wmin
            #println("Found $i ∈ I_wmin")
            f0 = fistar0(i0)
            #println("Compute fstar : $f")
            expo = count(j0 -> begin
                paddedInfoSet[j0+1] && inKi0(j0, i0, N0)
            end, f0:N0)
            #println("Computed exp $expo")
            c = enumarteSubtree(i0, i0, zeros(BIT, N))
            nweight += c*(1<<expo)
        end
    end
    return nweight
end

"""
    Note : this function should never be directly called, please prefer its alternative signature who saves compuation results
"""
function merge_rows(paddedInfoSet::Vector{Bool}, forbiden_set::Vector{NAT})
   N = length(paddedInfoSet)

    D = zeros(Bool, N)
    optD = zeros(Bool, N)

    R = Vector{NTuple{2,NAT}}([])
    optR = R
    
    Itrans = zeros(Bool, N)
    Fstar = zeros(Bool, N)
    Floser = zeros(Bool, N)
    Ftrans = Fstar


    dmin = minDist(paddedInfoSet, R)
    dmintrans= 0
    optDmin = 0

    admin = minWeightNumber(paddedInfoSet, R, dmin)
    admintrans = 0
    optAdmin = 0

    increaseWeight = true
    
    while true
        optDmin = dmin
        optAdmin = admin
        @LOG ROWMERGE_TAG "On main loop with increase_weight=$increaseWeight : opt_dmin=$optDmin, opdAdmin = $optAdmin"
        paddedWeightIndex!(N, optDmin, Itrans)
        for i ∈ eachindex(Itrans)
            if Itrans[i]
                increasingWeightFrozenRows!(paddedInfoSet, i, Fstar, Floser)
                all(j -> !j, Fstar) && continue

                if increaseWeight
                    Fstar .&= .!D
                    Ftrans = Fstar
                else
                    Floser .|= Fstar
                    Floser .&= .!D
                    Ftrans = Floser
                end
                c = 0
                for b in eachindex(Ftrans)
                    if c <= 2 && Ftrans[b] && b ∉ forbiden_set
                        b <= i && ERROR(ROWMERGE_TAG, "Should consider higher indices for merging rows")
                        c += 1
                        Rtrans = copy(R)
                        t::NTuple{2, NAT} = NTuple{2, NAT}((i,b))
                        push!(Rtrans, t)
                        @LOGD row_merge_tag "   Found merge candidate $t"
                        dmintrans = minDist(paddedInfoSet, Rtrans)
                        if dmintrans >= optDmin
                            @LOGD ROWMERGE_TAG "     Computing admin for a merge at index $b..."
                            admintrans = minWeightNumber(paddedInfoSet, Rtrans, dmintrans)
                            @LOG ROWMERGE_TAG "w=$dmintrans admin=$admintrans"
                            if dmintrans > optDmin || admintrans < optAdmin
                                @LOG ROWMERGE_TAG "Merges improve results : new opt_dmin=$dmintrans, new opt_admin=$admintrans"
                                optDmin = dmintrans
                                optAdmin = admintrans
                                copy!(optD, D)
                                optD[b] = true
                                optR = Rtrans
                            end
                        end
                    end
                end
            end
        end
        if optDmin > dmin || optAdmin < admin
            copy!(D,optD)
            R = optR
            dmin = optDmin
            admin = optAdmin
        elseif increaseWeight
            increaseWeight = false
        else
            @LOG ROWMERGE_TAG "Found $(length(R)) merges"
            return R
        end
    end
end
function merge_rows(infoset::Vector{NAT}, N::NAT, forbiden_set::Vector{NAT})
    id =  merges_already_present(N, infoset, forbiden_set)
    if id > 0 
        @LOG ROWMERGE_TAG "merges already computed"
        return merges_fetch(id)
    end
    padded_infoset = zeros(Bool, N)
    @inbounds for i in infoset
        padded_infoset[i] = true
    end
    merges = merge_rows(padded_infoset, forbiden_set)
    dump_merges(N, infoset, forbiden_set, merges)
    return merges
end
merge_rows(infoset::Vector{NAT}, N::NAT) = merge_rows(infoset, N, zeros(NAT, 0))

function rm_K(r::NAT, n::NAT)
    acc = 0
    for l ∈ 0:r
        acc += binomial(n,l)
    end
    return acc
end
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
function target_next_rm_code(N::NAT, K::NAT, reliabilities::Vector{REAL}, forbiden_set::Vector{NAT})
    ranking = sort(reverse(1:N), by=i->reliabilities[i])
    n = convert(NAT, log2(N))
    k = K
    while k <= N
        r = find_least_r(n, k)
        r < 0 && break 
        threshold = (1 << (n-r))
        rmInfoSet = filter(i -> polarRowWeigth0(i-1) >= threshold && i ∉ forbiden_set, ranking)
        length(rmInfoSet) >= K && return NAT.(rmInfoSet[1:K])
        k += 1
    end
    ERROR(ROWMERGE_TAG, "Couldn't find any merges for N=$N, K=$K")
end
target_next_rm_code(N::NAT, K::NAT, reliabilites::Vector{REAL}) = target_next_rm_code(N, K, reliabilites, zeros(NAT, 0))
