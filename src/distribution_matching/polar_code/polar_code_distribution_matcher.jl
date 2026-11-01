include("rr_sc.jl")

struct ConcretePolarCodeDistributionMatcher <: PolarDistributionMatcher
    _pc::ArbitraryPolarCode
    _use_rr::Bool
    _rr_scs::Vector{RandomizedRoundingSC}
    _init_llrs::Vector{REAL}
end

match!(d::ConcretePolarCodeDistributionMatcher, u::Vector{BIT}, dst::Vector{BIT}) = begin 

    d._use_rr ? polar_transform!(randomized_rounding_encode!(
        d._rr_scs[default_tid()],
        d._init_llrs,
        u,
        d._pc._padded_information_set,
        d._pc._padded_dynamically_shaping_set,
        d._pc._dynamically_shaping_set_func,
        dst
    )) : 
    encode!(d._pc, u, dst)
end
invert!(d::ConcretePolarCodeDistributionMatcher, x::Vector{BIT}, u::Vector{BIT}) = begin
    polar_transform!(x)
    cnt = 1
    @inbounds for i ∈ eachindex(x)
        if d._pc._padded_information_set[i]
            u[cnt] = x[i]
            cnt += 1
        end
    end
    return u
end
target_distribution(d::ConcretePolarCodeDistributionMatcher) = target_dist(d._pc)
info_rate(d::ConcretePolarCodeDistributionMatcher) = rate(d._pc)
blocklength(d::ConcretePolarCodeDistributionMatcher) = blocklength(d._pc)
messagelength(d::ConcretePolarCodeDistributionMatcher) = messagelength(d._pc)

function get_polar_code_distribution_matcher(args)
    blocklength = args["blocklength"]
    ispow2(blocklength) || ERROR(PCDM_TAG, "blocklength should be a power a two ($blocklength)")
    source_p = args["p"]
    0 <= source_p <= 1 || ERROR(PCDM_TAG, "$source_p is not a valid probability distribution")
    use_rr = args["randomized_rounding"]
    if use_rr
        args["encoding_list"] == 1 || ERROR(PCDM_TAG, "Cannot use randomized rounding when encoding list size is not 1 (it is $(args["encoding_list"]))")
    end
    hs = get_entropies_for_source(blocklength, source_p)
    source_index_reliabilites = sort(reverse(1:blocklength), by=i->hs[i])

    D = blocklength - floor(NAT, blocklength*args["shaping_rate"])

    if args["dynamic_shaping_constraint"] == "none"
        shaping_set = source_index_reliabilites[1:D]::Vector{NAT}
        length(shaping_set) == D || ERROR(PCDM_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
        dynamic_shaping_set = zeros(NAT, 0)
        shap_const = NoneConstraint()
    elseif args["dynamic_shaping_constraint"] == "crc"
        shaping_set = source_index_reliabilites[1:D]::Vector{NAT}
        length(shaping_set) == D || ERROR(PCDM_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
        sort!(shaping_set)
        dynamic_shaping_set = shaping_set[end-args["shaping_crc_width"]+1:end]::Vector{NAT}
        length(dynamic_shaping_set) == args["shaping_crc_width"] || ERROR(PCDM_TAG, "bug in construction (dynamic shaping set length should be $(args["shaping_crc_width"])), is $(length(dynamic_shaping_set))")
        filter!(i -> i ∉ dynamic_shaping_set, shaping_set)
        shap_const = get_crc_constraint(
            args["shaping_crc_width"],
            args["shaping_poly"],
            D-args["shaping_crc_width"],
            dynamic_shaping_set,
            args["shaping_known_at_dec"] ? ones(Bool, args["shaping_crc_width"]) : zeros(Bool, args["shaping_crc_width"])
        )
        length(shaping_set) + length(dynamic_shaping_set) == D || ERROR(PCDM_TAG, "bug in construction (shaping set length should be $(D-length(dynamic_shaping_set)), is $(length(shaping_set))")
    elseif args["dynamic_shaping_constraint"] == "rowmerge"
        #Warning : the selected shaping rate will not be the one acutally used
        shaping_set = target_next_rm_code(blocklength, D, hs)::Vector{NAT}
        length(shaping_set) == D || ERROR(PCDM_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
        merges = merge_rows(shaping_set, blocklength)
        dynamic_shaping_set = zeros(NAT, 0)
        for m in merges
            append!(dynamic_shaping_set, NAT(m[2]))
        end
        D += length(merges)
        shap_const = get_row_merge_constraint(merges)
    elseif args["dynamic_shaping_constraint"] == "paritycheck"
        k_check_points = args["shaping_k_pc_constraints"]
        shap_const, shaping_set = get_parity_check_contraints(
            k_check_points,
            hs,
            D/blocklength,
            min(args["encoding_list"]),
            true
        )
        dynamic_shaping_set = copy(shap_const._check_points)
        D += length(dynamic_shaping_set)
        length(dynamic_shaping_set) + length(shaping_set) == D || ERROR(PCDM_TAG, "bug in construction (shaping set length should be $(D-length(dynamic_shaping_set)), is $(length(shaping_set))")
    else
        ERROR(PCDM_TAG, "Unknown dynamic shaping bits constraint : $(args["dynamic_frozen_constraint"])")
    end

    info_set = filter(i -> i ∉ dynamic_shaping_set && i ∉ shaping_set, 1:blocklength)
    frozen_set = zeros(NAT, 0)
    dynamic_frozen_set = zeros(NAT, 0)
    froz_const = NoneConstraint()

    @LOGD PCDM_TAG "info_set_length : $(length(info_set)), dynamic_frozen_set_length : $(length(dynamic_frozen_set)), frozenset_length : $(length(frozen_set)), shaping_set_length : $(length(shaping_set)), dynamic_shaping_set_length : $(length(dynamic_shaping_set))"

    for set ∈ [shaping_set, dynamic_shaping_set, info_set, frozen_set, dynamic_frozen_set]
        for set2 ∈ [shaping_set, dynamic_shaping_set, info_set, frozen_set, dynamic_frozen_set]
            set != set2 && any(i -> i ∈ set2, set) && ERROR(PCDM_TAG, "sets with length $(length(set)), $(length(set2)) contain doublons")
        end
    end

    pc =  create_polar_code(
        blocklength,
        frozen_set,
        dynamic_frozen_set,
        shaping_set,
        dynamic_shaping_set,
        info_set,
        zeros(BIT, 0),
        args["encoding_list"],
        1,
        source_p,
        NoiselessChannel{BIT}(),
        froz_const,
        shap_const
    )
    init_llrs = zeros(REAL, blocklength)
    fill!(init_llrs, pc._source_llr)
    rr_sc = Vector{RandomizedRoundingSC}(undef, NDEFAULT)
    if use_rr
        @inbounds for i in 1:NDEFAULT
            rr_sc[i] = get_randomized_rounding_sc(blocklength)
        end
    end

    return ConcretePolarCodeDistributionMatcher(pc, use_rr, rr_sc, init_llrs)
end

Base.show(io::IO, d::DistributionMatcher) = begin 
    pc = d._pc
    pc._version_number == 1 ? 
        print(io, "PolarCode(blocklength=$(pc._blocklength);" *
                "shaping_set=$(NAT.(pc._padded_shaping_set));" *
                "dynamically_shaping_set=$(NAT.(pc._padded_dynamically_shaping_set));" *
                "information_set=$(NAT.(pc._padded_information_set));" *
                "encoding_list=$(pc._encoding_L);" *
                "source_p=$(target_dist(pc));" *
                "dynamic_shaping_bit_f=$(pc._dynamically_shaping_set_func);" *
                "randomized_rounding=$(NAT(d._use_rr))" *
                ")"
            ) : 
        print(io, "PolarCode(blocklength=$(pc._blocklength);" *
                "shaping_set=$(NAT.(pc._padded_shaping_set));" *
                "dynamically_shaping_set=$(NAT.(pc._padded_dynamically_shaping_set));" *
                "information_set=$(NAT.(pc._padded_information_set));" *
                "encoding_list=$(pc._encoding_L);" *
                "source_p=$(target_dist(pc));" *
                "dynamic_shaping_bit_f=$(pc._dynamically_shaping_set_func);" *
                "randomized_rounding=$(NAT(d._use_rr));" *
                "v=$(pc._version_number))"
        )
end
Base.parse(::Type{T}, s::AbstractString) where {T <: ConcretePolarCodeDistributionMatcher} = ABSTRACT(DM_TAG, T, "Base.parse")