struct ParityCheckConstraint <: DynamicConstraint
    _check_points::Vector{NAT}
    _constraints::Vector{Vector{NAT}}
end
get_parity_check_contraints(check_points::Vector{NAT}, constraints::Vector{Vector{NAT}}) = begin
    @inbounds for i in eachindex(check_points)
        k = check_points[i]
        @inbounds for v in constraints[i]
            v >= k && ERROR(PCDC_CONSTRAINT_PC_TAG, "partity check value $v is after computed check $k")
            v <= 0 && ERROR(PCDC_CONSTRAINT_PC_TAG, "partity check value $v is <= 0")
        end
    end
    return ParityCheckConstraint(check_points, constraints)
end

get_parity_check_contraints(k_check_points::NAT, hs::Vector{REAL}, rate::REAL, L::NAT, additional::Bool, forbiden_set::Vector{NAT},random::Bool) = begin
    0 <= rate <= 1 || ERROR(PCDC_CONSTRAINT_PC_TAG, "$rate is not within (0, 1)") 
    ispow2(L) || ERROR(PCDC_CONSTRAINT_PC_TAG, "$L should be a power of two")
    N = length(hs)
    ispow2(N) || ERROR(PCDC_CONSTRAINT_PC_TAG, "$N shoud be a poeer of two")
    l = max(1,convert(NAT, log2(L))-1)
    K = floor(NAT, rate*N)
    index_order = sort(reverse(1:N), by=i->hs[i])
    filter!(i -> i ∉ forbiden_set, index_order)
    if !random
        selected = additional ? index_order[K+1:K+k_check_points] : index_order[K-k_check_points+1:K]
        unselected = additional ? index_order[1:K] : index_order[1:K-k_check_points]
        length(selected) == k_check_points || ERROR(PCDC_CONSTRAINT_PC_TAG, "Internal bug")
        constraints = Vector{Vector{NAT}}(undef, k_check_points)
        for i in eachindex(constraints)
            k = selected[i]
            constraint = zeros(NAT, 0)
            for i in 1:l
                k-i >= 1 && append!(constraint, k - i)
            end
            constraints[i] = constraint
        end
        return get_parity_check_contraints(selected, constraints), unselected
    else
        selected = rand(index_order, k_check_points)
        filter!(i -> i ∉ selected, index_order)
        constraints = Vector{Vector{NAT}}(undef, k_check_points)
        for i in eachindex(constraints)
            k = selected[i]
            constraint = zeros(NAT, 0)
            for i in 1:l
                k-i >= 1 && append!(constraint, k - i)
            end
            constraints[i] = constraint
        end
        return get_parity_check_contraints(selected, constraints), additional ? index_order[1:K] : index_order[1:K-k_check_points]
    end
end

get_parity_check_contraints(k_check_points::NAT, hs::Vector{REAL}, rate::REAL, L::NAT, additional::Bool,random::Bool=true) = get_parity_check_contraints(k_check_points, hs, rate, L, additional, zeros(NAT, 0),random)


function get_check_point(P::ParityCheckConstraint, index::NAT)
    @inbounds for c in eachindex(P._check_points)
        P._check_points[c] == index && return c
    end
    ERROR(PCDC_CONSTRAINT_PC_TAG, "Called with wrong index")
end
(P::ParityCheckConstraint)(u::Vector{BIT}, select_index::Vector{Bool}, index::NAT) = begin
    index = get_check_point(P, index)
    acc = 0
    @inbounds for constraint_index in P._constraints[index]
        acc += u[constraint_index]
    end
    return BIT(acc % 2) 
end
known_at_decoding(P::ParityCheckConstraint, index::NAT) = true
Base.show(io::IO, P::ParityCheckConstraint) = print(io, "PCconstraint(check_points=$(P._check_points);constraints=$(P._constraints);randv=1)")
Base.parse(::Type{T}, s::AbstractString) where {T <: ParityCheckConstraint} = ABSTRACT(PCDC_CONSTRAINT_PC_TAG, T, "Base.parse")
clear_cache(P::ParityCheckConstraint) = nothing


