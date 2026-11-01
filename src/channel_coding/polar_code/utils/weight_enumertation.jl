function convert_i_to_bit_vector(i::NAT, dst::Vector{BIT})
           j = 1
           while i != 0
               dst[j] = (i&1)
               j += 1
               i >>= 1
           end
           return dst
       end

function brute_force_weigth_enum(N::NAT, K::NAT,dic)
    message = zeros(BIT, K)
    codeword = zeros(BIT, N)
    t = Dict{NAT, Vector{Vector{NAT}}}()
    for L in [1, 32]
        dic["encoding_list"] = L
        code = get_polar_code(dic)
        F = count(code._decoding_frozen_set)
        fv = zeros(BIT, F)
        weight_spectrum_for_fv = Vector{Vector{NAT}}(undef, (1 << F))
        weight_spectrum = zeros(NAT, N+1)
        @LOGD PC_UTILS_WEIGHT_ENUMERATOR "Starting weight enum for all $(1 << F) frozen vectors"
        for fv_i in 0:((1 << F)-1)
            fill!(weight_spectrum, NAT(0))
            #copying frozen_vec
            convert_i_to_bit_vector(fv_i, fv)
            for tid in 1:NDEFAULT
                fvc = code._statically_frozen_vec[tid]
                k = 1
                for i in eachindex(fvc)
                    if code._decoding_frozen_set[i]
                        fvc[i] = fv[k]
                        k += 1
                    end
                end
            end
            #exaustiv weight computation 
            @LOGD PC_UTILS_WEIGHT_ENUMERATOR "Starting weight enum for all $(1 << K) codewords with $fv_i frozen vector"
            for m_i in 0:((1 << K)-1)
                convert_i_to_bit_vector(m_i, message)
                encode!(code, message, codeword)
                cw_w = count(i->i==1, codeword)
                weight_spectrum[cw_w+1] += 1
            end
            weight_spectrum_for_fv[fv_i+1] = copy(weight_spectrum)
        end
        t[L] = weight_spectrum_for_fv
    end
    return t
end

function compute_avg_weight_for_frozen_vectors(N::NAT, weight_spectrum_for_fvs::Vector{Vector{NAT}})
    dic = zeros(REAL, length(weight_spectrum_for_fvs))
    multiplicator = [i for i in 0:N]
    for i in eachindex(weight_spectrum_for_fvs)
        w_s = weight_spectrum_for_fvs[i]
        avg_w = sum(w_s .* multiplicator)/(N*sum(w_s))
        dic[i] = avg_w
    end
    return dic
end

function diag_weight_transform(N, K)
    sigma = snr_db_to_sigma(DESIGN_SNR)
    optp,_ = opt_ook(sigma)
    W = OOK(N, sigma, optp)
    shap_rate = floor(NAT, N*(h2(optp) - 1/16)) / N
    dynam_shap_const = "paritycheck"
    dynam_froze_const = "none"

    dic = Dict(
        "blocklength" => N,
        "rate" => K/N,
        "design_channel" => W,
        "p" => optp, 
        "shaping_rate" => shap_rate,
        "dynamic_frozen_constraint" => dynam_froze_const,
        "frozen_crc_width" => 16,
        "frozen_poly" => get_default_poly(16),
        "frozen_known_at_dec" => true,
        "dynamic_shaping_constraint" => dynam_shap_const,
        "decoding_list" => 1,
        "encoding_list" => 1,
        "shaping_k_pc_constraints" => 10,
        "uniform_forzen_vec" => true
    )
    code = get_polar_code(dic)
    dic["encoding_list"] = 32
    code = get_polar_code(dic)
    t = brute_force_weigth_enum(N, K, dic)
    println("$(blocklength(code)), $(messagelength(code)), $(shaping_rate(code)), $(target_dist(code)), $(code._dynamically_frozen_set_func), $(code._dynamically_shaping_set_func)")
    for k in keys(t)
        avg1 = compute_avg_weight_for_frozen_vectors(N, t[k])
        avg_w2 = sum(avg1) / length(avg1)
        println("   L=$k ($avg_w2)")
        vec = t[k]
        fv = zeros(BIT, convert(NAT,log2(length(avg1))))
        for i in eachindex(vec)
            println("       $(avg1[i]) : $(vec[i]) $(convert_i_to_bit_vector(i-1,fv))")
        end
    end
end

