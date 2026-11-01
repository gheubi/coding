include("channel_code.jl")
include("simulation/channel_sim_data.jl")

function channel_simulate_once(code::Code{BIT,BIT,Y}, W::Channel{BIT, Y}, it::NAT) where {Y}
    it > 0 || ERROR(CHANNEL_SIM_TAG, "$it should be (stricly) greater than 0")  

    BEC = 0
    FEC = 0
    WEIGHT = 0
    WEIGHT2 = 0
    POW = 0.0
    NOISE = 0.0

    K = messagelength(code)
    N = blocklength(code)
    message = zeros(BIT, K)
    codeword = zeros(BIT, N)
    observation = zeros(Y, N)
    output_message = zeros(BIT, K)

    for i in 1:it
        @LOGD CHANNEL_SIM_TAG "============HERE WE GO AGAIN================"

        gen_uniform_binary_vec!(message)
        @LOGD CHANNEL_SIM_TAG "input message $(print_vector(message))"

        encode!(code, message, codeword)
        @LOGD CHANNEL_SIM_TAG "input codeword $(print_vector(codeword))"

        transmit!(W, codeword, observation)
        @LOGD CHANNEL_SIM_TAG "channel output $(print_vector(observation))"

        decode!(code, observation, output_message)
        @LOGD CHANNEL_SIM_TAG "decoder output $(print_vector(output_message))"

        #vec = code._decoding_full_vecs[default_tid()]
        #polar_transform!(vec)
        #WEIGHT2 += count(i -> i == 1, codeword)

        WEIGHT += count(i -> i == 1, codeword)
        bec = count(i -> output_message[i] != message[i], eachindex(message))
        bec > 0 && @LOGD CHANNEL_SIM_TAG "error at decoding : bec = $bec"
        BEC += bec
        FEC += (bec == 0) ? 0 : 1
        POW += last_signal_pow(W)
        NOISE += last_gen_noise_pow(W)
    end

    return (bec=BEC, fec=FEC, weight=WEIGHT, pow=POW, noise=NOISE, iters=it, w2=WEIGHT2)
end

function channel_simulate_par(code::Code{BIT, BIT, Y}, W::Channel{BIT,Y}, min_iters::NAT, max_iters::NAT, batch_size::NAT, already_computed_fec::NAT=0) where {Y}
    0 < min_iters <= max_iters || ERROR(CHANNEL_SIM_TAG, "Invalid relation between iterations numbers (min:$min_weight, max:$max_iters)")
    0 < batch_size || ERROR(CHANNEL_SIM_TAG, "$batch_size is not a valid batch_size")

    BEC = 0
    FEC = 0
    WEIGHT = 0
    #WEIGHT2 = 0
    POW = 0.0
    NOISE = 0.0

    frames = 0
    results = Vector{NamedTuple}(undef, NDEFAULT)


    while (already_computed_fec + FEC < 150) && (frames <= max_iters)
        @LOGD CHANNEL_SIM_TAG "About to call other threads"
        Threads.@threads for t ∈ 1:NDEFAULT
            results[t] = @time channel_simulate_once(code, W, batch_size)
        end
        @LOGD CHANNEL_SIM_TAG "All threads have returned"
        for tuple in results
            BEC += tuple.bec
            FEC += tuple.fec
            WEIGHT += tuple.weight
            POW += tuple.pow
            NOISE += tuple.noise
            frames += tuple.iters
            #WEIGHT2 += tuple.w2
        end
    end

    N = blocklength(code)
    AVG_WEIGTH = WEIGHT / (frames * N)
    AVG_SIG_POW = POW / (frames * N)
    AVG_NOISE_POW = NOISE / (frames * N)
    #AVG_WEIGTH2 = WEIGHT2 / (frames * N)
    #@LOG(CHANNEL_SIM_TAG, "encoding avg_weight=$AVG_WEIGTH, decoding avg_weight=$AVG_WEIGTH2")

    return (
        bec=BEC,
        fec=FEC,
        avg_weigth=AVG_WEIGTH,
        avg_sig_pow=AVG_SIG_POW,
        avg_noise_pow=AVG_NOISE_POW,
        iters=frames
    )
end

function sim_bawgn(
    code::Code{BIT, BIT, REAL},
    c::AbstractString,
    snrs_db::Vector{REAL},
    min_iters::NAT,
    max_iters::NAT,
    batch_size::NAT,
    change_chann::Bool,
    save_results::Bool=true
)
    c in ["BPSK", "OOK"] || ERROR(CHANNEL_SIM_TAG, "$c is an unknown binary constellation over the AWGNC")
    0 < min_iters <= max_iters || ERROR(CHANNEL_SIM_TAG, "Invalid relation between iterations numbers (min:$min_weight, max:$max_iters)")
    0 < batch_size || ERROR(CHANNEL_SIM_TAG, "$batch_size is not a valid batch_size")

    for snr_db in snrs_db
        @LOG CHANNEL_SIM_TAG "Simulating for blocklength=$(blocklength(code)) at snr=$snr_db 🤪"
        snr = value_from_db(snr_db)
        W = c == "OOK" ? OOK(blocklength(code), sqrt(1/snr), target_dist(code)) : BPSK(blocklength(code), sqrt(1/snr))
        change_chann && set_channel(code, W)
        fec = channel_code_get_computed_fec(code, c, change_chann, snr_db)
        r = channel_simulate_par(code, W, min_iters, max_iters, batch_size, fec)
        if save_results && r.iters > 0
            channel_code_dump_bin_results(
                code,
                c,
                snr_db,
                r.bec,
                r.fec,
                r.avg_weigth,
                r.avg_sig_pow,
                r.avg_noise_pow,
                r.iters,
                change_chann
            )
        else
            @LOG CHANNEL_SIM_TAG "at snr=$snr_db : r=$(rate(code)) sr=$(shaping_rate(code)) fec=$(r.fec) bec=$(r.bec) fer=$(r.fec/r.iters) avg_w=$(r.avg_weigth)"
        end
        @LOGD CHANNEL_SIM_TAG "done 😆"
    end

end