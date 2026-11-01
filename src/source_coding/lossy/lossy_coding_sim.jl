include("losssy_code.jl")

function lossy_coding_simulate_once(code::LossyCode{BIT,BIT,Y}, it::NAT) where {Y}
    it > 0 || ERROR(LOSSY_SIM_TAG, "$it should be (strictly) greater than 0")

    D = 0
    
    input_message = zeros(BIT, blocklength(code))
    encoded_message = zeros(BIT, messagelength(code))
    decoded_message = zeros(BIT, blocklength(code))
    p1 = input_distribution(code)
    for _ in 1:it
        gen_binary_vec!(p1, input_message)
        encode!(code, input_message, encoded_message)
        decode!(code, encoded_message, decoded_message)
        D += distortion(code, input_message, decoded_message)
    end
    return (d=D, frames=it)
end

include("simulation/lossy_code_fetch.jl")

function lossy_coding_simulate_par(code::LossyCode{BIT,BIT,Y}, target_dist::REAL, max_iters::NAT, batch_size::NAT) where {Y}
    0 < max_iters || ERROR(LOSSY_SIM_TAG, "$max_iters should be (stricly) greater than 0")
    0 < batch_size || ERROR(LOSSY_SIM_TAG, "$max_iters should be (stricly) greater than 0")

    @LOG LOSSY_SIM_TAG "Let's compute average distortion for rate $(rate(code)) (target $(target_dist)) 🤪"

    D = 0
    frames = 0
    results = Vector{NamedTuple}(undef, NDEFAULT)
    while frames <= max_iters
        Threads.@threads for t ∈ 1:NDEFAULT
            results[t] = @time lossy_coding_simulate_once(code, batch_size)
        end
        for tuple in results
            frames += tuple.frames
            D += tuple.d
        end
    end

    AVG_D = D / (frames * blocklength(code))
    dump_lossy_sim_res(code, target_dist, AVG_D, frames)

    @LOG LOSSY_SIM_TAG "Found average distortion $(AVG_D), it is $(AVG_D - target_dist) away from target 😆"
end
