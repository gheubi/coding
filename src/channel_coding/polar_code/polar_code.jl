include("polar_transform.jl")
include("scl.jl")
include("constructor/rowmerge.jl")
include("dynamic_constraint/dynamic_constraint.jl")

struct ArbitraryPolarCode{Y} <: PolarCode{Y}
    _blocklength::NAT
    _info_length::NAT

    _statically_frozen_vec::Vector{Vector{BIT}}

    _padded_dynamically_frozen_set::Vector{Bool}
    _padded_information_set::Vector{Bool}
    _padded_shaping_set::Vector{Bool}
    _padded_dynamically_shaping_set::Vector{Bool}

    _dynamically_frozen_set_func::DynamicConstraint
    _dynamically_shaping_set_func::DynamicConstraint

    _encoding_L::NAT
    _decoding_L::NAT

    _source_llr::REAL
    _design_channel::Channel{BIT,Y}
    _current_channel::Vector{Channel{BIT,Y}}

    _encoding_frozen_set::Vector{Bool}
    _encoding_dynamic_frozen_set::Vector{Bool}
    _decoding_frozen_set::Vector{Bool}
    _decoding_dynamic_frozen_set::Vector{Bool}

    _init_llrs::Vector{Vector{REAL}}
    _encoding_frozen_vecs::Vector{Vector{BIT}}
    _decoding_full_vecs::Vector{Vector{BIT}}
    _encoding_scl_working_ds::Vector{SCLWorkingDS}
    _decoding_scl_working_ds::Vector{SCLWorkingDS}
    _encoding_scl_res_ds::Vector{SCLResDS}
    _decoding_scl_res_ds::Vector{SCLResDS}
    _encoding_scl::Vector{SCLRes}
    _decoding_scl::Vector{SCLRes}

    _version_number::NAT
end

function decide_dynmically_frozen_bit(pc::ArbitraryPolarCode, u::Vector{BIT}, i::NAT)
    if pc._padded_dynamically_frozen_set[i]
        @LOGV POLAR_CODE_TAG "Dynamically choosing frozen bit $i..."
        return pc._dynamically_frozen_set_func(u, pc._padded_information_set, i)
    elseif pc._padded_dynamically_shaping_set[i]
        @LOGV POLAR_CODE_TAG "Dynamically choosing shaping bit $i..."
        return pc._dynamically_shaping_set_func(u, pc._padded_shaping_set, i)
    else
        ERROR(POLAR_CODE_TAG, "no dynamical decision to do at index $i")
    end
end
encode!(pc::ArbitraryPolarCode, u::Vector{BIT}, dst::Vector{BIT}) = begin
    if pc._version_number == 1
        encode_v1!(pc, u, dst)
    elseif pc._version_number == 2
        encode_v2!(pc, u, dst)
    else
        ERROR(POLAR_CODE_TAG, "Unknown version $(pc._version_number)")
    end
end
function encode_v2!(pc::ArbitraryPolarCode, u::Vector{BIT}, dst::Vector{BIT})
    @LOGV POLAR_CODE_TAG "encode!"
    @check length(dst) == pc._blocklength || ERROR(POLAR_CODE_TAG, "dst length's ($(length(dst))) does not match $(pc._blocklength)")
    @check length(u) == pc._info_length || ERROR(POLAR_CODE_TAG, "input message length ($(length(u))) does not match information set length $(pc._info_length)")
    encode_v1!(pc, u, dst)

    tid = default_tid()

    u_counter = 1
    embeded_vector = pc._encoding_frozen_vecs[tid]
    fill!(embeded_vector, 0)
    @inbounds for i in eachindex(dst)
        if pc._padded_information_set[i]
            embeded_vector[i] = u[u_counter]
            u_counter += 1
        elseif pc._decoding_frozen_set[i]
            embeded_vector[i] = pc._statically_frozen_vec[tid][i]
        elseif pc._padded_dynamically_frozen_set[i] #&& !known_at_decoding(pc._dynamically_frozen_set_func, i)
            embeded_vector[i] = pc._dynamically_frozen_set_func(embeded_vector, pc._padded_information_set, i)
        elseif pc._padded_dynamically_shaping_set[i]
            embeded_vector[i] = pc._dynamically_shaping_set_func(embeded_vector, pc._padded_shaping_set, i)
        end
    end
    is_valid_message(pc, embeded_vector) || ERROR(POLAR_CODE_TAG, "bug in encoding")
    @LOGD POLAR_CODE_TAG "embeded input $(print_vector(embeded_vector))"
    llrs = pc._init_llrs[tid]
    @inbounds for i in eachindex(llrs)
        llrs[i] = pc._source_llr + LLR(pc._current_channel[1], transmit_one_symbol(pc._current_channel[1], dst[i]))
    end

    scl = scl_decode!(
        llrs,
        pc._encoding_L,
        pc,
        pc._encoding_scl_working_ds[tid],
        pc._encoding_scl_res_ds[tid],
        pc._encoding_frozen_set,
        embeded_vector,
        pc._encoding_dynamic_frozen_set,
        pc._encoding_scl[tid]
    )
    mostlikely!(scl, pc, dst)
    return polar_transform!(dst)
end
function encode_v1!(pc::ArbitraryPolarCode, u::Vector{BIT}, dst::Vector{BIT})
    @LOGV POLAR_CODE_TAG "encode!"
    @check length(dst) == pc._blocklength || ERROR(POLAR_CODE_TAG, "dst length's ($(length(dst))) does not match $(pc._blocklength)")
    @check length(u) == pc._info_length || ERROR(POLAR_CODE_TAG, "input message length ($(length(u))) does not match information set length $(pc._info_length)")


    tid = default_tid()

    u_counter = 1
    embeded_vector = pc._encoding_frozen_vecs[tid]
    fill!(embeded_vector, 0)
    @inbounds for i in eachindex(dst)
        if pc._padded_information_set[i]
            embeded_vector[i] = u[u_counter]
            u_counter += 1
        elseif pc._decoding_frozen_set[i]
            embeded_vector[i] = pc._statically_frozen_vec[tid][i]
        elseif pc._padded_dynamically_frozen_set[i] #&& !known_at_decoding(pc._dynamically_frozen_set_func, i)
            embeded_vector[i] = pc._dynamically_frozen_set_func(embeded_vector, pc._padded_information_set, i)
        elseif pc._padded_dynamically_shaping_set[i]
            embeded_vector[i] = pc._dynamically_shaping_set_func(embeded_vector, pc._padded_shaping_set, i)
        end
    end
    is_valid_message(pc, embeded_vector) || ERROR(POLAR_CODE_TAG, "bug in encoding")
    @LOGD POLAR_CODE_TAG "embeded input $(print_vector(embeded_vector))"
    fill!(pc._init_llrs[tid], pc._source_llr)
    if shaping_rate(pc) != 1.0
        @LOGV POLAR_CODE_TAG "calling scl to encode..."
        scl = scl_decode!(
            pc._init_llrs[tid],
            pc._encoding_L,
            pc,
            pc._encoding_scl_working_ds[tid],
            pc._encoding_scl_res_ds[tid],
            pc._encoding_frozen_set,
            embeded_vector,
            pc._encoding_dynamic_frozen_set,
            pc._encoding_scl[tid]
        )
        @LOGD POLAR_CODE_TAG "scl-encoded $scl"
        mostlikely!(scl, pc, dst)
        @LOGV POLAR_CODE_TAG "found most likely encoded vector"
    else
        @LOGV POLAR_CODE_TAG "no need to use scl to encode"
        copy!(dst, embeded_vector)
    end
    return polar_transform!(dst)
end

function decode!(pc::ArbitraryPolarCode{Y}, y::Vector{Y}, dst::Vector{BIT}) where {Y}
    @LOGV POLAR_CODE_TAG "decode!"
    @check length(y) == pc._blocklength || ERROR(POLAR_CODE_TAG, "input length's ($(length(y))) does not match $(pc._blocklength)")
    @check length(dst) == pc._info_length || ERROR(POLAR_CODE_TAG, "message length ($(length(dst))) does not match information set length $(pc._info_length)")

    tid = default_tid()

    llrs = pc._init_llrs[tid]
    @inbounds for i in eachindex(y)
        llrs[i] = pc._source_llr + LLR(pc._current_channel[1],y[i]) 
    end
    @LOGD POLAR_CODE_TAG "Init llrs : $(print_vector(llrs))"
    @LOGV POLAR_CODE_TAG "calling scl to decode..."
    @LOGD POLAR_CODE_TAG "before decoding : scl metrics are $(print_vector(pc._decoding_scl[tid].metrics))"
    scl = scl_decode!(
        llrs,
        pc._decoding_L,
        pc,
        pc._decoding_scl_working_ds[tid],
        pc._decoding_scl_res_ds[tid],
        pc._decoding_frozen_set,
        pc._statically_frozen_vec[tid],
        pc._decoding_dynamic_frozen_set,
        pc._decoding_scl[tid]
    )
    @LOGD POLAR_CODE_TAG "after decoding : scl metrics are $(print_vector(pc._decoding_scl[tid].metrics))"
    res = mostlikely!(scl, pc, pc._decoding_full_vecs[tid])
    @LOGD POLAR_CODE_TAG "scl-decoded \n$scl"
    @LOGD POLAR_CODE_TAG "best message : $(print_vector(res))"

    @LOGV POLAR_CODE_TAG "found most likely decoded vector"
    dst_counter = 1
    @inbounds for i in eachindex(res)
        if pc._padded_information_set[i]
            dst[dst_counter] = res[i]
            dst_counter += 1
        end
    end
    @LOGV POLAR_CODE_TAG "compressed most likely vector into message"
    clear_cache(pc._dynamically_frozen_set_func)
    clear_cache(pc._dynamically_shaping_set_func)
    return dst
end

blocklength(c::ArbitraryPolarCode) = c._blocklength
messagelength(c::ArbitraryPolarCode) = c._info_length
set_channel(c::ArbitraryPolarCode{Y}, W::Channel{BIT,Y}) where {Y} = begin
    c._current_channel[1] = W
end
is_valid_message(c::ArbitraryPolarCode, u::Vector{BIT}) = begin
    @LOGV POLAR_CODE_TAG "is_valid_message have been called "
    length(u) == c._blocklength || return false
    tid = default_tid()
    @inbounds for i in eachindex(u)
        if c._decoding_frozen_set[i]
            u[i] == c._statically_frozen_vec[tid][i] || return false
        elseif c._padded_dynamically_frozen_set[i]
            u[i] == c._dynamically_frozen_set_func(u, c._padded_information_set, i) || return false 
        elseif c._padded_dynamically_shaping_set[i]
            u[i] == c._dynamically_shaping_set_func(u, c._padded_shaping_set, i) || return false
        else
            (u[i] == 0 || u[i] == 1) || return false
        end
    end
    return true
end
is_valid_codeword(c::ArbitraryPolarCode, x::Vector{BIT}) = begin
    polar_transform!(x)
    res = is_valid_message(c, x)
    polar_transform!(x)
    return res
end
rate(c::ArbitraryPolarCode) = c._info_length / c._blocklength
shaping_rate(c::ArbitraryPolarCode) = (c._blocklength - count(c._padded_dynamically_shaping_set) - count(c._padded_shaping_set)) / c._blocklength
target_dist(c::ArbitraryPolarCode) = 1/(exp(c._source_llr)+1)
current_channel(c::ArbitraryPolarCode) = c._current_channel[1]

##########################################################################################################
###### Some utilitaries function  ########################################################################
##########################################################################################################

include("constructor/polar_code_constuctor.jl")
include("bit_channel_entropies/find_polarized_entropies.jl")
include("utils/weight_enumertation.jl")

get_indices_back(paddeed_v::Vector{Bool}) = begin
    set = zeros(NAT, 0)
    for i in eachindex(paddeed_v)
        if paddeed_v[i]
            append!(set, i)
        end
    end
    return NAT.(set)
end

print_readable_pc(pc::ArbitraryPolarCode) = begin
    info_set = get_indices_back(pc._padded_information_set)
    frozen_set = get_indices_back(pc._decoding_frozen_set)
    dyn_frozen_set = get_indices_back(pc._padded_dynamically_frozen_set)
    shap_set = get_indices_back(pc._padded_shaping_set)
    dyn_shap_set = get_indices_back(pc._padded_dynamically_shaping_set)
    println("infoset : $info_set")
    println("frozenset : $frozen_set")
    println("dynamicly_frozen : $dyn_frozen_set")
    println("shapingset : $shap_set")
    println("dyn_shapingset : $dyn_shap_set")
    println("frozen_set_f : $(pc._dynamically_frozen_set_func)")
    println("shaping_set_f : $(pc._dynamically_shaping_set_func)")
end

to_tikz(pc::ArbitraryPolarCode) = begin
    info_set = get_indices_back(pc._padded_information_set)
    frozen_set = get_indices_back(pc._decoding_frozen_set)
    dyn_frozen_set = get_indices_back(pc._padded_dynamically_frozen_set)
    shap_set = get_indices_back(pc._padded_shaping_set)
    dyn_shap_set = get_indices_back(pc._padded_dynamically_shaping_set)
    s = ""
    for i in 1:pc._blocklength
        if i ∈ info_set
            color = "white!70!blue"
        elseif i ∈ frozen_set
            color = "white!80"
        elseif i ∈ dyn_frozen_set
            color = "white!70!yellow"
        elseif i ∈ shap_set
            color = "white!70!orange"
        elseif i ∈ dyn_shap_set
            color = "white!50!orange"
        else
            ERROR(POLAR_CODE_TAG, "unknown index $i")
        end
        s1 = "\\fill[$color] ($(i-1), 0) rectangle ($i, 0.8);\n"        
        s2 = "\\draw  ($(i-1), 0) rectangle ($i, 0.8);\n"
        s *= s1 
        s *= s2
    end
    return s
end


Base.show(io::IO, pc::ArbitraryPolarCode) = begin
    pc._version_number == 1 ? 
        print(io, "PolarCode(blocklength=$(pc._blocklength);" *
                "statically_frozen_set=$(NAT.(pc._decoding_frozen_set));" *
                "dynamically_frozen_set=$(NAT.(pc._padded_dynamically_frozen_set));" *
                "shaping_set=$(NAT.(pc._padded_shaping_set));" *
                "dynamically_shaping_set=$(NAT.(pc._padded_dynamically_shaping_set));" *
                "information_set=$(NAT.(pc._padded_information_set));" *
                "statically_frozen_vec=$(NAT.(pc._statically_frozen_vec[1]));" *
                "encoding_list=$(pc._encoding_L);" *
                "decoding_list=$(pc._decoding_L);" *
                "source_p=$(target_dist(pc));" *
                "design_channel=$(pc._design_channel);" *
                "dynamic_frozen_bit_f=$(pc._dynamically_frozen_set_func);" *
                "dynamic_shaping_bit_f=$(pc._dynamically_shaping_set_func)" *
                ")"
            ) : 
        print(io, "PolarCode(blocklength=$(pc._blocklength);" *
                "statically_frozen_set=$(NAT.(pc._decoding_frozen_set));" *
                "dynamically_frozen_set=$(NAT.(pc._padded_dynamically_frozen_set));" *
                "shaping_set=$(NAT.(pc._padded_shaping_set));" *
                "dynamically_shaping_set=$(NAT.(pc._padded_dynamically_shaping_set));" *
                "information_set=$(NAT.(pc._padded_information_set));" *
                "statically_frozen_vec=$(NAT.(pc._statically_frozen_vec[1]));" *
                "encoding_list=$(pc._encoding_L);" *
                "decoding_list=$(pc._decoding_L);" *
                "source_p=$(target_dist(pc));" *
                "design_channel=$(pc._design_channel);" *
                "dynamic_frozen_bit_f=$(pc._dynamically_frozen_set_func);" *
                "dynamic_shaping_bit_f=$(pc._dynamically_shaping_set_func)" *
                "v=$(pc._version_number))"
        )
end


Base.parse(::Type{T}, s::AbstractString) where {T <: PolarCode} = begin
    startswith(s, "PolarCode(") || ERROR(POLAR_CODE_TAG, "input string cannot be parsed")
    dic = parse_to_dic(
        s,
        POLAR_CODE_TAG,
        [
        "blocklength",
        "statically_frozen_set",
        "dynamically_frozen_set",
        "shaping_set",
        "dynamically_shaping_set",
        "information_set",
        "statically_frozen_vec",
        "encoding_list",
        "decoding_list",
        "source_p",
        "design_channel",
        "dynamic_frozen_bit_f",
        "dynamic_shaping_bit_f",
        ]
    )
    blocklength = parse(NAT, dic["blocklength"])
    statically_frozen_set = parse_vector(Bool, dic["statically_frozen_set"])
    dynamically_frozen_set = parse_vector(Bool, dic["dynamically_frozen_set"])
    shaping_set = parse_vector(Bool, dic["shaping_set"])
    dynamically_shaping_set = parse_vector(Bool, dic["dynamically_shaping_set"])
    information_set = parse_vector(Bool, dic["information_set"])
    statically_frozen_vec = parse_vector(BIT, dic["statically_frozen_vec"])
    encoding_list = parse(NAT, dic["encoding_list"])
    decoding_list = parse(NAT, dic["decoding_list"])
    source_p = parse(REAL, dic["source_p"])
    channel = parse(Channel, dic["design_channel"])
    dynamic_frozen_bit_f = parse(DynamicConstraint, dic["dynamic_frozen_bit_f"])
    dynamic_shaping_bit_f = parse(DynamicConstraint, dic["dynamic_shaping_bit_f"])
    version_number = "v" in keys(dic) ? parse(NAT, dic["v"]) : 1

    return create_polar_code(
        blocklength,
        statically_frozen_set,
        dynamically_frozen_set,
        shaping_set,
        dynamically_shaping_set,
        information_set,
        statically_frozen_vec,
        encoding_list,
        decoding_list,
        source_p,
        channel,
        dynamic_frozen_bit_f,
        dynamic_shaping_bit_f,
        version_number
    )
end