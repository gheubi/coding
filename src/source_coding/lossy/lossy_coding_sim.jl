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


#=
    Simulate a lossy code by encoding and decoding random input messages and computing the average distortion.
    The function takes a lossy code, a target distortion, a maximum number of iterations, and a batch size as input.
    It generates random input messages according to the input distribution of the code, encodes them, decodes them, and computes the distortion between the original and reconstructed messages.
    The simulation is performed in parallel using multiple threads to speed up the computation.
    The results are dumped to a CSV file for later analysis.
    @param code the lossy code to simulate
    @param target_dist the target distortion for the simulation
    @param max_iters the maximum number of iterations to run the simulation
    @param batch_size the number of messages to process in each iteration
    @return the average distortion, and the results are saved to a CSV file    
=#
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
    return AVG_D
end
