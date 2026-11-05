using CSV, DataFrames

# checks is a similar polar code has already been simulates
# two polar codes are similar if they agree up to the frozen bits values
# (except all zero frozen vector, which is considered to be different from any other frozen vector values)
function polar_code_already_present_up_to_fv(code::PolarCode{T}) where {T}
    c = "$code"
    c_chunks = split(c, ";statically_frozen_vec=[")
    c_chunkss = split(c_chunks[2],';')[2:end]

    df = CSV.read("data/simulation/channel_coding_params.csv", DataFrame)
    for row in eachrow(df)
        r_chunks = split(row.code, ";statically_frozen_vec=[")
        length(r_chunks) != 2 && continue
        r_chunkss = split(r_chunks[2], ';')[2:end]
        if c_chunks[1] == r_chunks[1] && c_chunkss == r_chunkss && row.buggy == 0
            @LOGD POLAR_CODE_CONSTRUCTOR_TAG "Found $(row.code) to be quite similar to $c"
            return (row.code, row.id)
        end
    end
    return ("", -1)
end

# extracts the frozen vector from a polarcode string
function get_frozen_vec(s::AbstractString)
    s = split(split(s, "statically_frozen_vec=")[2], ';')[1]
    fv = parse_vector(BIT, s)
    return fv
end

#=
    creates a polar code

    @param blocklength : the blocklength of the polar code, should be a power of 2
    @param statically_frozen_set : the set of statically frozen bits
    @param dynamically_frozen_set : the set of dynamically frozen bits
    @param shaping_set : the set of shaping bits
    @param dynamically_shaping_set : the set of dynamically shaping bits
    @param information_set : the set of information bits
    @param statically_frozen_vec : the value of the statically frozen bits
    @param encoding_list : the encoding list size, should be a power of 2
    @param decoding_list : the decoding list size, should be a power of 2
    @param source_p : the target binary-distribution of the codewords
    @param channel : the design channel used for the code constructiom
    @param dynamic_frozen_bit_f : the structure describing the computation of the dynamic frozen bits
    @param dynamic_shaping_bit_f : the structure describing the computation of the dynamic shaping bits
    @param version_number : the version number of the polar code, default 1 (1 is the only supported version for now)

    @return a polar code with the specified parameters
=#
function create_polar_code(
    blocklength::NAT,
    statically_frozen_set::Vector{Bool},
    dynamically_frozen_set::Vector{Bool},
    shaping_set::Vector{Bool},
    dynamically_shaping_set::Vector{Bool},
    information_set::Vector{Bool},
    statically_frozen_vec::Vector{BIT},
    encoding_list::NAT,
    decoding_list::NAT,
    source_p::REAL,
    channel::Channel{BIT,Y},
    dynamic_frozen_bit_f::DynamicConstraint,
    dynamic_shaping_bit_f::DynamicConstraint,
    version_number::NAT=1
) where {Y}
    @LOG POLAR_CODE_CONSTRUCTOR_TAG "Creating polar code with blocklength $blocklength"
    ispow2(blocklength) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Blocklength must be a power of two for polar codes ($blocklength)")
    0 <= source_p <= 1 || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "$source_p is not a valid probability")
    length(statically_frozen_set) == length(dynamically_frozen_set) == 
        length(shaping_set) == length(dynamically_shaping_set) == 
        length(information_set) == length(statically_frozen_vec) == blocklength ||
        ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set don't have a size equal to blocklength")
    (ispow2(encoding_list) && ispow2(decoding_list)) || 
        ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "encoding list size ($encoding_list) or decoding list size ($decoding_list) is not a power of two")

    all((NAT.(statically_frozen_set) .+ NAT.(dynamically_frozen_set) .+ NAT.(shaping_set) .+
        NAT.(dynamically_shaping_set) .+ NAT.(information_set)) .== 1) ||
        ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition")

    (ispow2(encoding_list) && ispow2(decoding_list)) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "encoding or decoding list size are not a power of two")
        
    
    source_llr = log(1-source_p) - log(source_p)

    encoding_frozen_set = zeros(Bool, blocklength)
    encoding_dynamic_frozen_set = zeros(Bool, blocklength)
    decoding_frozen_set = zeros(Bool, blocklength)
    decoding_dynamic_frozen_set = zeros(Bool, blocklength)
    @inbounds for i in eachindex(encoding_frozen_set)
        encoding_frozen_set[i] = information_set[i] || statically_frozen_set[i]
        encoding_dynamic_frozen_set[i] = 
            (dynamically_frozen_set[i] && known_at_decoding(dynamic_frozen_bit_f, i)) || 
            (dynamically_shaping_set[i] && known_at_decoding(dynamic_shaping_bit_f, i))
        decoding_frozen_set[i] = statically_frozen_set[i]
        decoding_dynamic_frozen_set[i] = 
            (dynamically_frozen_set[i] && known_at_decoding(dynamic_frozen_bit_f, i)) || 
            (dynamically_shaping_set[i] && known_at_decoding(dynamic_shaping_bit_f, i))
    end

    init_llrs = Vector{Vector{REAL}}(undef, NDEFAULT)
    encoding_frozen_vecs = Vector{Vector{BIT}}(undef, NDEFAULT)
    decoding_full_vecs = Vector{Vector{BIT}}(undef, NDEFAULT)
    encoding_scl_working_ds = Vector{SCLWorkingDS}(undef, NDEFAULT)
    decoding_scl_working_ds = Vector{SCLWorkingDS}(undef, NDEFAULT)
    encoding_scl_res_ds = Vector{SCLResDS}(undef, NDEFAULT)
    decoding_scl_res_ds = Vector{SCLResDS}(undef, NDEFAULT)
    encoding_scl = Vector{SCLRes}(undef, NDEFAULT)
    decoding_scl = Vector{SCLRes}(undef, NDEFAULT)
    statically_frozen_vecs = Vector{Vector{BIT}}(undef, NDEFAULT)
    @inbounds for i in 1:NDEFAULT
        init_llrs[i] = zeros(REAL, blocklength)
        encoding_frozen_vecs[i] = zeros(BIT, blocklength)
        decoding_full_vecs[i] = zeros(BIT, blocklength)
        encoding_scl_working_ds[i] = get_scl_working_ds(blocklength, encoding_list)
        encoding_scl_working_ds[i].ptr == C_NULL && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "ALlocation error from C")
        decoding_scl_working_ds[i] = get_scl_working_ds(blocklength, decoding_list)
        decoding_scl_working_ds[i].ptr == C_NULL && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "ALlocation error from C")
        encoding_scl_res_ds[i] = get_scl_res_ds(blocklength, encoding_list)
        encoding_scl_res_ds[i].ptr == C_NULL && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "ALlocation error from C")
        decoding_scl_res_ds[i] = get_scl_res_ds(blocklength, decoding_list)
        decoding_scl_res_ds[i].ptr == C_NULL && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "ALlocation error from C")
        encoding_scl[i] =  emptySCL(blocklength, encoding_list)
        decoding_scl[i] = emptySCL(blocklength, decoding_list)
        statically_frozen_vecs[i] = copy(statically_frozen_vec)
    end

    pc = ArbitraryPolarCode{Y}(
        blocklength,
        count(information_set),
        statically_frozen_vecs,
        dynamically_frozen_set,
        information_set,
        shaping_set,
        dynamically_shaping_set,
        dynamic_frozen_bit_f,
        dynamic_shaping_bit_f,
        encoding_list,
        decoding_list,
        source_llr,
        channel,
        [channel],
        encoding_frozen_set,
        encoding_dynamic_frozen_set,
        decoding_frozen_set,
        decoding_dynamic_frozen_set,
        init_llrs,
        encoding_frozen_vecs,
        decoding_full_vecs,
        encoding_scl_working_ds,
        decoding_scl_working_ds,
        encoding_scl_res_ds,
        decoding_scl_res_ds,
        encoding_scl,
        decoding_scl,
        version_number
    )
    s, pid = polar_code_already_present_up_to_fv(pc)
    pid == -1 && return pc
    @LOG POLAR_CODE_CONSTRUCTOR_TAG "fetching old polar code with id $pid"
    fv = get_frozen_vec(s)
    @inbounds for i in 1:NDEFAULT
        copy!(pc._statically_frozen_vec[i], fv)
    end
    return pc
end

#=
    creates a polar code

    @param blocklength : the blocklength of the polar code, should be a power of 2
    @param statically_frozen_set : the set of statically frozen bits
    @param dynamically_frozen_set : the set of dynamically frozen bits
    @param shaping_set : the set of shaping bits
    @param dynamically_shaping_set : the set of dynamically shaping bits
    @param information_set : the set of information bits
    @param statically_frozen_vec : the value of the statically frozen bits
    @param encoding_list : the encoding list size, should be a power of 2
    @param decoding_list : the decoding list size, should be a power of 2
    @param source_p : the target binary-distribution of the codewords
    @param channel : the design channel used for the code constructiom
    @param dynamic_frozen_bit_f : the structure describing the computation of the dynamic frozen bits
    @param dynamic_shaping_bit_f : the structure describing the computation of the dynamic shaping bits
    @param version_number : the version number of the polar code, default 1 (1 is the only supported version for now)

    @return a polar code with the specified parameters
=#
function create_polar_code(
    blocklength::NAT,
    statically_frozen_set::Vector{NAT},
    dynamically_frozen_set::Vector{NAT},
    shaping_set::Vector{NAT},
    dynamically_shaping_set::Vector{NAT},
    information_set::Vector{NAT},
    statically_frozen_vec::Vector{BIT},
    encoding_list::NAT,
    decoding_list::NAT,
    source_p::REAL,
    channel::Channel{BIT,Y},
    dynamic_frozen_bit_f::DynamicConstraint,
    dynamic_shaping_bit_f::DynamicConstraint,
    version_number::NAT=1
) where {Y}
    ispow2(blocklength) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Blocklength must be a power of two for polar codes ($blocklength)")
    0 <= source_p <= 1 || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "$source_p is not a valid probability")
    (ispow2(encoding_list) && ispow2(decoding_list)) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "encoding or decoding list size are not a power of two")

    length(unique(statically_frozen_set)) == length(statically_frozen_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set contains duplicates")
    length(unique(dynamically_frozen_set)) == length(dynamically_frozen_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set contains duplicates")
    length(unique(shaping_set)) == length(shaping_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set contains duplicates")
    length(unique(dynamically_shaping_set)) == length(dynamically_shaping_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set contains duplicates")
    length(unique(information_set)) == length(information_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input set contaidynamic_shaping_bit_fns duplicates")
    length(statically_frozen_vec) == length(statically_frozen_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input frozen vector has wrong size")

    sum = length(statically_frozen_set) + length(dynamically_frozen_set) +  length(shaping_set) + length(dynamically_shaping_set) + length(information_set) 
    sum == blocklength || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets have wrong size (total size is $sum : should be $blocklength)")
    
    all_in = zeros(Bool, blocklength)
    padded_statically_frozen_set = zeros(Bool, blocklength)
    padded_frozen_vec = zeros(BIT, blocklength)
    fc = 1
    @inbounds for i in statically_frozen_set
        padded_statically_frozen_set[i] = true
        padded_frozen_vec[i] = statically_frozen_vec[fc]
        fc += 1
        all_in[i] && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
        all_in[i] = true
    end
    padded_dynamically_frozen_set = zeros(Bool, blocklength)
    @inbounds for i in dynamically_frozen_set
        padded_dynamically_frozen_set[i] = true
        all_in[i] && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
        all_in[i] = true
    end
    padded_shaping_set = zeros(Bool, blocklength)
    @inbounds for i in shaping_set
        padded_shaping_set[i] = true
        all_in[i] && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
        all_in[i] = true
    end
    padded_dynamically_shaping_set = zeros(Bool, blocklength)
    @inbounds for i in dynamically_shaping_set
        padded_dynamically_shaping_set[i] = true
        all_in[i] && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
        all_in[i] = true
    end
    padded_infoset = zeros(Bool, blocklength)
    @inbounds for i in information_set
        padded_infoset[i] = true
        all_in[i] && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
        all_in[i] = true
    end
    all(all_in) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Input sets don't form a partition of {1, ..., blocklength}")
    
    return create_polar_code(
        blocklength,
        padded_statically_frozen_set,
        padded_dynamically_frozen_set,
        padded_shaping_set,
        padded_dynamically_shaping_set,
        padded_infoset,
        padded_frozen_vec,
        encoding_list,
        decoding_list,
        source_p,
        channel,
        dynamic_frozen_bit_f,
        dynamic_shaping_bit_f,
        version_number
    )
end

#=
    creates a polar code

    @param args : a dic containing the values to create the polar code, it should contain the following keys 
        - "blocklength" : the blocklength of the code, should be a power of 2
        - "design_channel" : the design channel used to construct the code by selecting the reliable bit-channels
        - "p" : the target bitwise binary-distribution of the codewords
        - "shaping_rate" : the shaping rate of the code (1 - #number of bits used for shaping/blocklength)
        - "dynamic_shaping_constraint" : the type of dynamic shaping constraint, should be among :
            - "none" : no dynamic shaping constraint
            - "crc" : use crc-check, if used, the following keys must be added:
                - "shaping_crc_width" : the number of crc bits
                - "shaping_poly" : the crc binary polynomial used to compute the remainder
                - "shaping_known_at_decoding" : determines whether the dynamic shaping bits are computed during the SCL
                    computation, or after in order to select the most-likely valid path
            - "rowmerge" : use row-merging. Uf used, the algortihm will select a set of merges and the shaping rate will be 
                adapted based on the number of merges found. In that case, the resulting code will have a different shaping
                rate as the one specified in args.
            - "paritycheck" : the k indices are selected to be encoded as the binary sum of the log2(min(encdoingL, decodingL)).
                preceeding bits. If used, the following key must be added:
                - "shaping_k_pc_constraints" : the number of parity check bits
                The indices for partiy check bits are found as the most-relibable (with respect to source polarization), non-used for shaping bits
        - "encoding_list" : the encoding list size, should be a power of 2
        - "rate" : the communication rate of the code
        - "dynamic_frozen_constraint" : the type if dynamic frozen bits constraint, should be among :
            - "none" : no dynamic frozen bits constraint
            - "crc" : use crc-check, if used, the following keys must be added:
                - "frozen_crc_width" : the number of crc bits
                - "frozen_poly" : the crc binary polynomial used to compute the remainder
                - "frozen_known_at_decoding" : determines whether the dynamic frozen bits are computed during the SCL
                    computation, or after in order to select the most-likely valid path
            - "rowmerge" : use row-merging. Uf used, the algortihm will select a set of merges and the rate will be 
                adapted based on the number of merges found. In that case, the resulting code will have a different rate
                rate as the one specified in args.
                @warning the rowmerge procedure will likely fail to found optimal merges if there are some shaping bits and/or
                    dynamical shaping bits.
            - "paritycheck" : the k indices are selected to be encoded as the binary sum of the log2(min(encdoingL, decodingL)).
                preceeding bits. If used, the following key must be added:
                - "info_k_pc_constraints" : the number of parity check bits
                The indices for partiy check bits are found as the most-relibable (with respect to channel polarization), non-used for information bits
        - "uniform_forzen_vec" : if false, the frozen vector is all 0, otherwise it is set uniformly at random.

    @return a polar code with the specified parameters
=#
function get_polar_code(args)
    blocklength = args["blocklength"]
    ispow2(blocklength) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "blocklength should be a power a two ($blocklength)")
    design_w = args["design_channel"]
    source_p = args["p"]
    0 <= source_p <= 1 || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "$source_p is not a valid probability distribution")

    hs, hw = get_entropies_for_source_and_channel(blocklength, source_p, design_w)
    source_index_reliabilites = sort(1:blocklength, by=i->hs[i])
    channel_index_reliabilites = sort(reverse(1:blocklength), by=i->hw[i])

    D = ceil(NAT, blocklength*(1-args["shaping_rate"]))
    if args["dynamic_shaping_constraint"] == "none"
        shaping_set = source_index_reliabilites[1:D]::Vector{NAT}
        length(shaping_set) == D || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
        dynamic_shaping_set = zeros(NAT, 0)
        shap_const = NoneConstraint()
    elseif args["dynamic_shaping_constraint"] == "crc"
        shaping_set = source_index_reliabilites[1:D]::Vector{NAT}
        length(shaping_set) == D || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
        sort!(shaping_set)
        dynamic_shaping_set = shaping_set[end-args["shaping_crc_width"]+1:end]::Vector{NAT}
        length(dynamic_shaping_set) == args["shaping_crc_width"] || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in construction (dynamic shaping set length should be $(args["shaping_crc_width"])), is $(length(dynamic_shaping_set))")
        filter!(i -> i ∉ dynamic_shaping_set, shaping_set)
        shap_const = get_crc_constraint(
            args["shaping_crc_width"],
            args["shaping_poly"],
            D-args["shaping_crc_width"],
            dynamic_shaping_set,
            args["shaping_known_at_dec"] ? ones(Bool, args["shaping_crc_width"]) : zeros(Bool, args["shaping_crc_width"])
        )
        length(shaping_set) + length(dynamic_shaping_set) == D || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in construction (shaping set length should be $(D-length(dynamic_shaping_set)), is $(length(shaping_set))")
    elseif args["dynamic_shaping_constraint"] == "rowmerge"
        #Warning : the selected shaping rate will not be the one acutally used
        shaping_set = target_next_rm_code(blocklength, D, hs)::Vector{NAT}
        length(shaping_set) == D || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (shaping set length should be $D, is $(length(shaping_set)))")
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
            min(args["encoding_list"], args["decoding_list"]),
            true
        )
        dynamic_shaping_set = copy(shap_const._check_points)
        D += length(dynamic_shaping_set)
        length(dynamic_shaping_set) + length(shaping_set) == D || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in construction (shaping set length should be $(D-length(dynamic_shaping_set)), is $(length(shaping_set))")
    else
        ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Unknown dynamic shaping bits constraint : $(args["dynamic_frozen_constraint"])")
    end

    K = floor(NAT,blocklength*args["rate"])
    if args["dynamic_frozen_constraint"] == "none"
        F = blocklength - (K+D)
        frozen_set = channel_index_reliabilites[end-F+1:end]::Vector{NAT}
        info_set = filter(i -> i ∉ shaping_set && i ∉ frozen_set && i ∉ dynamic_shaping_set, 1:blocklength)
        dynamic_frozen_set = zeros(NAT, 0)
        info_set = filter(i -> i ∉ shaping_set && i ∉ dynamic_shaping_set && i ∉ frozen_set && i ∉ dynamic_frozen_set, 1:blocklength)::Vector{NAT}
        length(info_set) == K || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (info set length should be $K, is $(length(info_set)))")
        froz_const = NoneConstraint()
    elseif args["dynamic_frozen_constraint"] == "crc"
        K += args["frozen_crc_width"]
        F = blocklength-(K+D)
        frozen_set = channel_index_reliabilites[end-F+1:end]::Vector{NAT}
        length(frozen_set) == F || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (frozen set length should be $F, is $(length(frozen_set))")
        info_set = filter(i -> i ∉ shaping_set && i ∉ frozen_set && i ∉ dynamic_shaping_set, 1:blocklength)
        length(info_set) == K || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (info set length should be $K, is $(length(info_set)))")
        sort!(info_set)
        dynamic_frozen_set = info_set[end-args["frozen_crc_width"]+1:end]::Vector{NAT}
        length(dynamic_frozen_set) == args["frozen_crc_width"] || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (dynmic frozen set length should be $(args["frozen_crc_width"]), is $(length(dynamic_frozen_set)))")
        filter!(i -> i ∉ dynamic_frozen_set, info_set)
        froz_const = get_crc_constraint(
            args["frozen_crc_width"],
            args["frozen_poly"],
            K - args["frozen_crc_width"],
            dynamic_frozen_set,
            args["frozen_known_at_dec"] ? ones(Bool, args["frozen_crc_width"]) : zeros(Bool, args["frozen_crc_width"])
        )
    elseif args["dynamic_frozen_constraint"] == "rowmerge"
        info_set = target_next_rm_code(blocklength, K, hw, [shaping_set; dynamic_shaping_set])::Vector{NAT}
        length(info_set) == K || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (info set length should be $K, is $(length(info_set)))")
        merges = merge_rows(info_set, blocklength, [shaping_set; dynamic_shaping_set])
        dynamic_frozen_set = zeros(NAT,0)
        for m in merges
            append!(dynamic_frozen_set, NAT(m[2]))
        end
        frozen_set = filter(i -> i ∉ shaping_set && i ∉ dynamic_shaping_set && i ∉ info_set && i ∉ dynamic_frozen_set, 1:blocklength)::Vector{NAT}
        froz_const = get_row_merge_constraint(merges)
    elseif args["dynamic_frozen_constraint"] == "paritycheck"
        k_check_points = args["info_k_pc_constraints"]
        froz_const, info_set = get_parity_check_contraints(
            k_check_points,
            hw,
            K/blocklength,
            args["decoding_list"],
            true,
            [dynamic_shaping_set; shaping_set]
        )
        dynamic_frozen_set = copy(froz_const._check_points)
        length(info_set) == K || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "bug in constructing (info set length should be $K, is $(length(info_set)))")
        frozen_set = filter(i -> i ∉ shaping_set && i ∉ dynamic_shaping_set && i ∉ info_set && i ∉ dynamic_frozen_set, 1:blocklength)::Vector{NAT}
        @check all(i -> 1 <= i <= blocklength, frozen_set) || ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "frozen_set contains unvalid indices : $(filter(i -> i <= 0 || i >= blocklength, frozen_set))")
        @check any(i -> i ∈ shaping_set, frozen_set) && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "frozen set contains indices in shaping set")
        @check any(i -> i ∈ dynamic_shaping_set, frozen_set) && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "frozen set contains indices in dynamic shaping set")
        @check any(i -> i ∈ info_set, frozen_set) && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "frozen set contains indices in info set")
        @check any(i -> i ∈ dynamic_frozen_set, frozen_set) && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "frozen set contains indices in dynmamic frozen set")
    else
        ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "Unknown dynamic frozen bits constraint : $(args["dynamic_frozen_constraint"])")
    end

    @LOGD POLAR_CODE_CONSTRUCTOR_TAG "info_set_length : $(length(info_set)), dynamic_frozen_set_length : $(length(dynamic_frozen_set)), frozenset_length : $(length(frozen_set)), shaping_set_length : $(length(shaping_set)), dynamic_shaping_set_length : $(length(dynamic_shaping_set))"

    for set ∈ [shaping_set, dynamic_shaping_set, info_set, frozen_set, dynamic_frozen_set]
        for set2 ∈ [shaping_set, dynamic_shaping_set, info_set, frozen_set, dynamic_frozen_set]
            set != set2 && any(i -> i ∈ set2, set) && ERROR(POLAR_CODE_CONSTRUCTOR_TAG, "sets with length $(length(set)), $(length(set2)) contain doublons")
        end
    end

    return create_polar_code(
        blocklength,
        NAT.(frozen_set),
        NAT.(dynamic_frozen_set),
        NAT.(shaping_set),
        NAT.(dynamic_shaping_set),
        NAT.(info_set),
        args["uniform_forzen_vec"] ? gen_uniform_binary_vec(BIT, length(frozen_set)) : zeros(BIT, length(frozen_set)),
        args["encoding_list"],
        args["decoding_list"],
        source_p,
        design_w,
        froz_const,
        shap_const
    )
end