include("fetch_polarized_entropies.jl")


function index_entropies_simulate_once(pc::PolarCode{Y}, p1::REAL, W::Channel{BIT, Y}, it::NAT) where {Y}
    tid = default_tid()
    N = blocklength(pc)
    u = zeros(BIT, messagelength(pc))
    y = zeros(Y, N)
    dumy = zeros(BIT, blocklength(pc))

    entropies = zeros(REAL, N)

    for _ in 1:it
        @LOGV POLARIZE_ENTROPY_TAG "Generating binary vector"
        gen_binary_vec!(p1, dumy)

        @LOGV POLARIZE_ENTROPY_TAG "Simulating transmission"
        transmit!(W, dumy, y)

        polar_transform!(dumy)

        begin #this part shouldn't be interupted btw

            tid = default_tid() #updating tid (the damn scheduler can always be messing with us 👀)
            #genie aided decoding
            copy!(pc._statically_frozen_vec[tid], dumy)
            @LOGV POLARIZE_ENTROPY_TAG "decoding..."
            decode!(pc, y, u)

            @LOGV POLARIZE_ENTROPY_TAG "Adding llrs"
            llrs = pc._decoding_scl[tid].llr[1]

        end #now can be interupted 

        @inbounds for i in eachindex(llrs)
            #println( h2(1/(1+exp(llrs[i]))))
            entropies[i] += h2(1/(1+exp(llrs[i])))
        end
    end

    return (hs=entropies, frames=it)
end

function index_entropies_simulate_par(blocklength::NAT, p1::REAL, W::Channel{BIT, Y}, iters::NAT, batch_size::NAT) where {Y}
    ispow2(blocklength) || ERROR(POLARIZE_ENTROPY_TAG, "$blocklength should be a power of 2")
    0 <= p1 <= 1 || ERROR(POLARIZE_ENTROPY_TAG, "$p1 is not a valid probability")

    pc = create_polar_code(
        blocklength,
        ones(Bool, blocklength),
        zeros(Bool, blocklength),
        zeros(Bool, blocklength),
        zeros(Bool, blocklength),
        zeros(Bool, blocklength),
        zeros(BIT, blocklength),
        1,
        1,
        p1,
        W,
        NoneConstraint(),
        NoneConstraint()
    )

    @LOG POLARIZE_ENTROPY_TAG "Constructing fake polar code : done"

    res = Vector{NamedTuple}(undef, NDEFAULT)
    frames = 0
    h = zeros(REAL, blocklength)
    while frames < iters
        Threads.@threads for t ∈ 1:NDEFAULT
            res[t] = @time index_entropies_simulate_once(pc, p1, W, batch_size)
        end
        for tuple in res
            @inbounds for i in eachindex(h)
                h[i] += tuple.hs[i]
            end
            frames += tuple.frames 
        end
    end

    h ./= frames

    return (h, frames)
end

function index_entropies_simulate_source(blocklength::NAT, p1::REAL, iters::NAT, batch_size::NAT)
    @LOG POLARIZE_ENTROPY_TAG "Computing polarized indexes for source with $p1..."
    nl = NoiselessChannel{BIT}()
    hs, its = index_entropies_simulate_par(blocklength, p1, nl, iters, batch_size)
    entropy_dump_results(blocklength, p1, nl, hs, its)
    @LOG POLARIZE_ENTROPY_TAG "done"
    return hs
end

function index_entropies_simulate_channel(blocklength::NAT, p1::REAL, W::Channel{BIT, Y}, iters::NAT, batch_size::NAT) where {Y}
    @LOG POLARIZE_ENTROPY_TAG "Computing polarized indexes for source with $p1 and channel $W..."
    hw, itw = index_entropies_simulate_par(blocklength, p1, W, iters, batch_size)
    entropy_dump_results(blocklength, p1, W, hw, itw)
    @LOG POLARIZE_ENTROPY_TAG "done"
    return hw
end

function index_entropies_simulate_source_and_chann(blocklength::NAT, p1::REAL, W::Channel{BIT, Y}, iters::NAT, batch_size::NAT) where {Y}
    
    hs = index_entropies_simulate_source(blocklength, p1, iters, batch_size)
    hw = index_entropies_simulate_channel(blocklength, p1, W, iters, batch_size)

    return hs, hw
end

function get_entropies_for_source(blocklength::NAT, p1::REAL)
    nl = NoiselessChannel{BIT}()
    hs = entropy_data_already_present(blocklength, p1, nl) > 0 ? 
            entropy_fetch_hs(blocklength, p1, nl) : index_entropies_simulate_source(blocklength, p1, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE)
    return hs
end

function get_entropies_for_channel(blocklength::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    hw = entropy_data_already_present(blocklength, p1, W) > 0 ? 
            entropy_fetch_hs(blocklength, p1, W) : index_entropies_simulate_channel(blocklength, p1, W, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE)
    return hw
end

function get_entropies_for_source_and_channel(blocklength::NAT, p1::REAL, W::Channel{BIT, Y}) where {Y}
    hs = get_entropies_for_source(blocklength, p1)
    hw = get_entropies_for_channel(blocklength, p1, W)
    return hs, hw
end

function gen_a_lot_of_entropies(Ns=[4,8,16,32,64,128,256,512,1024], p1=[i/8 for i in 1:7])
    sigma = snr_db_to_sigma(DESIGN_SNR)
    for N in Ns, p in p1
        W = OOK(N, sigma, p)
        index_entropies_simulate_source_and_chann(N, p, W, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE)
        W = BPSK(N, sigma)
        index_entropies_simulate_source_and_chann(N, p, W, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE)
    end
end

function gen_entropies_for_opt_p(Ns=[4,8,16,32,64,128,256,512,1024], snrs_db = [i/8 for i in -32:64])
    for snr_db in snrs_db
        sigma = snr_db_to_sigma(snr_db)
        opt_p, opt_i = opt_ook(sigma)
        for N in Ns
            W = OOK(N, sigma, opt_p)
            index_entropies_simulate_source_and_chann(N, opt_p, W, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE)
        end
    end
end