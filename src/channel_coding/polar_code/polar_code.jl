include("polar_transform.jl")
include("scl.jl")
include("constructor/rowmerge.jl")
include("dynamic_constraint/dynamic_constraint.jl")

struct ArbitraryPolarCode{Y} <: PolarCode{Y}
    _blocklength::NAT                                       # the blocklength of the code N
    _info_length::NAT                                       # the length of the information set K

    _statically_frozen_vec::Vector{Vector{BIT}}             # the padded vector of statically frozen bits (the frozen bits that are known at both encoding and decoding)

    _padded_dynamically_frozen_set::Vector{Bool}            # padded vector of dynamically frozen bits (indicates which indices are dynamically frozen) 
                                                            # a dynamically frozen bit is a bit that is frozen at encoding but not known at decoding (it is determined by a function of the (preceeding) information bits)
    _padded_information_set::Vector{Bool}                   # padded vector of information bits (indicates which indices are information bits)
    _padded_shaping_set::Vector{Bool}                       # padded vector of shaping bits (indicates which indices are shaping bits)
                                                            # a shaping bit is a bit used to target a specific distribution of the codeword (it is determined by a function of the (preceeding) bits)
    _padded_dynamically_shaping_set::Vector{Bool}           # padded vector of dynamically shaping bits (indicates which indices are used for dynamic shaping)
                                                            # a dynamically shaping bit is a bit that can be used both for shaping and forward error correction, its goal is to obtain a better shaping 
                                                            # performance than by only using the shaping set. The functions used for the dynamically shaping bits are typically linear and always different
                                                            # than the one used for the shaping bits, which relies on the most likely bits given the preceeding ones.

    _dynamically_frozen_set_func::DynamicConstraint         # dynamic frozen bit structure, used to determine the value of a dynamically frozen bit based on the preceeding bits
    _dynamically_shaping_set_func::DynamicConstraint        # dynamic shaping bit structure, used to determine the value of a dynamically shaping bit based on the preceeding bits

    _encoding_L::NAT                                        # the size of the list used for successive cancellation list (SCL) encoding (shaping)
    _decoding_L::NAT                                        # the sizeo of the list used for SCL decoding (shaping and fec)

    _source_llr::REAL                                       # the log-likelihood ratio (LLR) of the target source distribution (used for shaping)
    _design_channel::Channel{BIT,Y}                         # the channel used to design the code (select the reliable and unreliable bits)
    _current_channel::Vector{Channel{BIT,Y}}                # vector (size 1) of channel used to compute the probability of the received symbols (used for decoding), can be modified and/or different from the design channel 

    _encoding_frozen_set::Vector{Bool}                      # padded vector of the bits considered as frozen at encoding (statically frozen bit + information bits)
    _encoding_dynamic_frozen_set::Vector{Bool}              # padded vector of the bits considered as dynamic at encoding (dynamically frozen bits, dynamically shaping bits) : some bits are technically dynamic but not necessariky 
                                                            # computed during the SCL encoding (aka known_at_deocding). Instead, they are used to choose the final codeword among the candidates in the SCL list. 
    _decoding_frozen_set::Vector{Bool}                      # padded vector of the bits considered as frozen at decoding (statically frozen bit)  
    _decoding_dynamic_frozen_set::Vector{Bool}              # padded vector of the bits considered as dynamic at decoding (dynamically frozen bits, dynamically shaping bits) : some bits are technically dynamic but not necessariky 
                                                            # computed during the SCL decoding (aka known_at_deocding). Instead, they are used to choose the final codeword among the candidates in the SCL list. 

    _init_llrs::Vector{Vector{REAL}}                        # vector (size of the number of threads) of vectors storing the LLRs of the bit-channels, used both at encoding and decoding
    _encoding_frozen_vecs::Vector{Vector{BIT}}              # vector (size of the number of threads) of vectors storing the bits considered as frozen at encoding (statically frozen bit + information bits)
    _decoding_full_vecs::Vector{Vector{BIT}}                # vector (size of the number of threads) of vectors used to store the most-likely and valid codeword found after SCL both at encoding and decoding
    _encoding_scl_working_ds::Vector{SCLWorkingDS}          # vector (size of the number of threads) of working data structures used for SCL encoding
    _decoding_scl_working_ds::Vector{SCLWorkingDS}          # vector (size of the number of threads) of working data structures used for SCL decoding (not the same DS as for encoding, as the list size can be different)
    _encoding_scl_res_ds::Vector{SCLResDS}                  # vector (size of the number of threads) of result data structures used for SCL encoding    
    _decoding_scl_res_ds::Vector{SCLResDS}                  # vector (size of the number of threads) of result data structures used for SCL decoding (not the same DS as for encoding, as the list size can be different)
    _encoding_scl::Vector{SCLRes}                           # vector (size of the number of threads) of SCL structure result used for SCL encoding
    _decoding_scl::Vector{SCLRes}                           # vector (size of the number of threads) of SCL structure result used for SCL decoding (not the same DS as for encoding, as the list size can be different)

    _version_number::NAT                                    # version number of the polar code implementation, used to select the appropriate encoding and decoding functions and keep compatibility as well as history in the simulation results
end

#=
    Determine the value of a dynamically frozen bit based on the preceeding bits and the information set.
    The function takes an ArbitraryPolarCode, a vector of bits u, and an index i as input.
    It checks if the bit at index i is dynamically frozen or dynamically shaping, and calls the appropriate function to determine its value.
    If the bit is not dynamically frozen or shaping, it raises an error.
    @param pc the ArbitraryPolarCode containing the dynamic frozen and shaping set functions
    @param u the vector of bits used to determine the value of the dynamically frozen or shaping bit
    @param i the index of the bit to determine
    @return the value of the dynamically frozen or shaping bit at index i
=#
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
#=
    Encode a vector of information bits u into a vector of output bits dst.
    The function takes an ArbitraryPolarCode, a vector of information bits u, and a vector of output bits dst as input.
    It checks that the lengths of u and dst match the information set length and blocklength of the code, respectively.
    It then calls the appropriate encoding function based on the version number of the code.
    @param pc the ArbitraryPolarCode to use for encoding
    @param u the vector of information bits to encode
    @param dst the vector to store the encoded output bits
    @return dst
=#
encode!(pc::ArbitraryPolarCode, u::Vector{BIT}, dst::Vector{BIT}) = begin
    if pc._version_number == 1
        return encode_v1!(pc, u, dst)
    elseif pc._version_number == 2
        return encode_v2!(pc, u, dst)
    else
        ERROR(POLAR_CODE_TAG, "Unknown version $(pc._version_number)")
    end
end
#=
    Encode a vector of information bits u into a vector of output bits dst using version 2 of the encoding function.
    This version is deprecated due to poor performance and should not be used.
    The function takes an ArbitraryPolarCode, a vector of information bits u, and a vector of output bits dst as input.
    It checks that the lengths of u and dst match the information set length and blocklength of the code, respectively.
    It then embeds the information bits into the appropriate positions in the codeword, computes the log-likelihood ratios (LLRs) for each bit-channel, 
    and performs successive cancellation list (SCL) encoding (with dynamic frozen and shaping bits) to find the most likely codeword.
    Finally, it applies the polar transform to the output bits.
    @param pc the ArbitraryPolarCode to use for encoding
    @param u the vector of information bits to encode
    @param dst the vector to store the encoded output bits
    @return dst
=#
function encode_v2!(pc::ArbitraryPolarCode, u::Vector{BIT}, dst::Vector{BIT})
    ERROR(POLAR_CODE_TAG, "encode_v2! is deprecated, please use version 1 instead : version 2 suffers from terrible performance")
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
        #those two cases are not necessary, as the indices will be recomputed
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
        #test computing the probabilites as if the message to transmmit had already been encoded and transmited
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
#=
    Encode a vector of information bits u into a vector of output bits dst using version 1 of the encoding function.
    The function takes an ArbitraryPolarCode, a vector of information bits u, and a vector of output bits dst as input.
    It checks that the lengths of u and dst match the information set length and blocklength of the code, respectively.
    It then embeds the information bits into the appropriate positions in the codeword, computes the log-likelihood ratios (LLRs) for each bit-channel, 
    and performs successive cancellation list (SCL) encoding (with dynamic frozen and shaping bits) to find the most likely codeword.
    Finally, it applies the polar transform to the output bits.
    @param pc the ArbitraryPolarCode to use for encoding
    @param u the vector of information bits to encode
    @param dst the vector to store the encoded output bits
    @return dst
=#
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

#=
    Decode a vector of received symbols y into a vector of reconstructed bits dst.
    The function takes an ArbitraryPolarCode, a vector of received symbols y, and a vector of reconstructed bits dst as input.
    It checks that the lengths of y and dst match the blocklength and information set length of the code, respectively.
    It then computes the log-likelihood ratios (LLRs) for each bit-channel based on the received symbols and the current channel,
    and performs successive cancellation list (SCL) decoding (with dynamic frozen and shaping bits) to find the most likely codeword.
    Finally, it extracts the information bits from the decoded codeword and stores them in dst.
    @param pc the ArbitraryPolarCode to use for decoding
    @param y the vector of received symbols to decode
    @param dst the vector to store the reconstructed information bits
    @return dst
=#
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

# returns the blocklength of the code
blocklength(c::ArbitraryPolarCode) = c._blocklength

# returns the messagelength of the code
messagelength(c::ArbitraryPolarCode) = c._info_length

# set the channel used to compute the probability of the received symbols (used for decoding), can be modified and/or different from the design channel
set_channel(c::ArbitraryPolarCode{Y}, W::Channel{BIT,Y}) where {Y} = begin
    c._current_channel[1] = W
end

#=
    Check if a given vector of bits u is a valid embedded message for the given ArbitraryPolarCode.
    The function takes an ArbitraryPolarCode and a vector of bits u as input.
    It checks that the length of u matches the blocklength of the code, and then iterates through each bit in u to verify that it satisfies the constraints imposed by the frozen sets and dynamic frozen/shaping sets.
    If any bit does not satisfy the constraints, the function returns false. Otherwise, it returns true.
    @param c the ArbitraryPolarCode to check against
    @param u the vector of bits to check
    @return true if u is a valid embedded message, false otherwise
=#
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
#=
    Check if a given vector of bits x is a valid codeword for the given ArbitraryPolarCode.
    The function takes an ArbitraryPolarCode and a vector of bits x as input.
    It applies the polar transform to x, checks if the resulting vector is a valid embedded message using is_valid_message, and then applies the polar transform again to restore x to its original state.
    If the transformed vector is a valid message, the function returns true. Otherwise, it returns false.
    @param c the ArbitraryPolarCode to check against
    @param x the vector of bits to check
    @return true if x is a valid codeword, false otherwise
=#
is_valid_codeword(c::ArbitraryPolarCode, x::Vector{BIT}) = begin
    polar_transform!(x)
    res = is_valid_message(c, x)
    polar_transform!(x)
    return res
end

# returns the rate of the polar code
rate(c::ArbitraryPolarCode) = c._info_length / c._blocklength

# returns the shaping rate of the polar code, which is the fraction of bits not used for shaping 1 - (shaping bits + dynamically shaping bits) / blocklength
shaping_rate(c::ArbitraryPolarCode) = (c._blocklength - count(c._padded_dynamically_shaping_set) - count(c._padded_shaping_set)) / c._blocklength

# returns the target distribution of the polar code, which is the probability of a bit being 1 in the codeword
target_dist(c::ArbitraryPolarCode) = 1/(exp(c._source_llr)+1)

# returns the channel used to compute the probability of the received symbols (used for decoding), can be modified and/or different from the design channel
current_channel(c::ArbitraryPolarCode) = c._current_channel[1]

##########################################################################################################
###### Some utilitaries function  ########################################################################
##########################################################################################################

include("constructor/polar_code_constuctor.jl")
include("bit_channel_entropies/find_polarized_entropies.jl")
include("utils/weight_enumertation.jl")

# tranforms a padded vector of booleans into a vector of indices (1-based) where the boolean is true
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

#= 
    Generate a TikZ representation of the ArbitraryPolarCode, which can be used to visualize the code structure in LaTeX documents.
    The function takes an ArbitraryPolarCode as input and returns a string containing the TikZ code.
    It identifies the indices of the information set, frozen set, dynamically frozen set, shaping set, and dynamically shaping set, and assigns colors to each type of bit.
    The resulting TikZ code consists of filled rectangles representing each bit in the codeword, with different colors indicating their roles in the code.
    @param pc the ArbitraryPolarCode to visualize
    @return a string containing the TikZ code for the visualization
=#
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
    ERROR(POLAR_CODE_TAG, "parsing PolarCode is not implemented yet")
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