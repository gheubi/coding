include("distribution_matcher.jl")
include("simulation/distribution_matcher_fetch.jl")
function distribution_matching_simulate_once(dm::DistributionMatcher{BIT}, it::NAT) 
    it > 0 || ERROR(DM_SIM_TAG, "$it should be (strictly) greater than 0")

    p1 = 0
    errors = 0

    u = zeros(BIT, messagelength(dm))
    x = zeros(BIT, blocklength(dm))
    û = zeros(BIT, messagelength(dm))

    for _ in 1:it
        gen_uniform_binary_vec!(u)
        match!(dm, u, x)
        #=
        polar_transform!(x)
        cnt = 1
        @inbounds for i in eachindex(x)
            if dm._pc._padded_information_set[i]
                u[cnt] = x[i]
                cnt += 1
            end
        end
        match!(dm, u, x)
        =#
        p1 += weight(x)
        invert!(dm,x,û)
        if u != û
            errors += 1
        end
    end

    return (ONES=p1, ERRORS=errors, FRAMES=it)
end

function distribution_matching_simulate_par(dm::DistributionMatcher{BIT}, max_iters::NAT, batch_size::NAT,dump::Bool=true)
    0 < max_iters || ERROR(DM_SIM_TAG, "$max_iters should be (stricly) greater than 0")
    0 < batch_size || ERROR(DM_SIM_TAG, "$max_iters should be (stricly) greater than 0")

    @LOG DM_SIM_TAG "Let's compute average distribution for info rate $(info_rate(dm)) (target $(target_distribution(dm))) 🤪"

    ones = 0
    fec = 0
    frames = 0
    results = Vector{NamedTuple}(undef, NDEFAULT)
    while frames <= max_iters
        Threads.@threads for t ∈ 1:NDEFAULT
            results[t] = @time distribution_matching_simulate_once(dm, batch_size)
        end
        for tuple in results
            frames += tuple.FRAMES
            ones += tuple.ONES
            fec += tuple.ERRORS
        end
    end

    avg_one = ones / (frames * blocklength(dm))
    fer = fec / frames

    if dump
        dump_distribution_matching_sim_res(dm, avg_one, fer, frames)
    else
        @LOG DM_SIM_TAG "Found $avg_one, $fer"
    end

    @LOG DM_SIM_TAG "Found averge distribution $(avg_one) with fer $(fer), it is $(avg_one - target_distribution(dm)) away from target 😆"
end