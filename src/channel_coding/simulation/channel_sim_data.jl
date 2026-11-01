using CSV, DataFrames

const channel_params_path = "src/channel_coding/simulation/data/channel_coding_params.csv"
const channel_sim_path = "src/channel_coding/simulation/data/channel_coding_sim.csv"
const channel_dm_path = "src/channel_coding/simulation/data/distribution_matching_sim.csv"

function channel_code_already_present(code::Code{U,X,REAL}, constellation::AbstractString, change_chann::Bool) where {U,X}
    constellation in ["OOK", "BPSK"] || ERROR(CHANNEL_SIM_DATA_TAG, "$constellation is not a known constellation")

    df = CSV.read(channel_params_path, DataFrame)
    cs = "$code"
    bl = blocklength(code)
    r = rate(code)
    sr = shaping_rate(code)
    cc = Int(change_chann)
    for row in eachrow(df)
        row.blocklength == bl && row.rate == r && row.shaping_rate == sr && 
            row.code == cs && row.constellation == constellation && row.buggy == 0 && row.channel_changed == cc &&
            return row.id
    end
    return -1
end

function channel_code_dump_params(code::Code{U,X,REAL}, constellation::AbstractString, change_chann::Bool) where {U,X}
    constellation in ["OOK", "BPSK"] || ERROR(CHANNEL_SIM_DATA_TAG, "$constellation is not a known constellation")
    @LOGD CHANNEL_SIM_DATA_TAG "Checking whether polar code is already in"
    id = channel_code_already_present(code, constellation, change_chann::Bool)
    id > 0 && return id
    @LOGD CHANNEL_SIM_DATA_TAG "It isn't"
    df = CSV.read(channel_params_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, blocklength(code), rate(code), shaping_rate(code), constellation, "$code", 0, Int(change_chann)], promote=true)
    CSV.write(channel_params_path, df)
    return id
end

function channel_code_get_computed_fec(code::Code{U,X,REAL}, constellation::AbstractString, change_chann::Bool, snr::REAL) where {U,X}
    id = channel_code_already_present(code, constellation, change_chann)
    id <= 0 && return NAT(0)
    df = CSV.read(channel_sim_path, DataFrame)
    for row in eachrow(df)
        if row.id == id && row.snr == snr 
            @LOG CHANNEL_SIM_DATA_TAG "For $id at snr $snr : already computed fec $(row.fec)"
            return NAT(row.fec)
        end
    end
    return NAT(0)
end

function channel_code_params_id_alread_in(id::NAT, df::DataFrame)
    any(df.id .== id)
end

function channel_code_dump_sim_res(id::NAT, snr::REAL, bec::NAT, fec::NAT, frames::NAT, message_len::NAT)

    ber = bec / (message_len * frames)
    fer = fec / frames

    0 <= fer <= 1 || ERROR(CHANNEL_SIM_DATA_TAG, "$fer is not in the validity range [0, 1]")
    0 <= ber <= 1 || ERROR(CHANNEL_SIM_DATA_TAG, "$ber is not in the validity range [0, 1]")



    df = CSV.read(channel_sim_path, DataFrame)
    found = false
    if channel_code_params_id_alread_in(id, df)
        for row in eachrow(df)
            if row.id == id && row.snr == snr
                row.bec += bec
                row.fec += fec
                row.frames += frames
                row.ber = row.bec / (message_len * row.frames)
                row.fer = row.fec / row.frames
                found = true
                break
            end
        end
    end
    if !found
        push!(df, [id, snr, fer, ber, frames, fec, bec], promote=true)
    end
    CSV.write(channel_sim_path, df)
    return id
end

function channel_code_dump_bin_shap_res(id::NAT, shap_rate::REAL, target_dist::REAL, target_snr::REAL, observed_dist::REAL, avg_sig_pow::REAL, avg_noise_pow::REAL, frames::NAT)
    df = CSV.read(channel_dm_path, DataFrame)
    found = false
    if channel_code_params_id_alread_in(id, df)
        for row in eachrow(df)
            if row.id == id && row.target_snr == target_snr && row.target_dist == target_dist && row.shaping_rate == shap_rate
                new_frames = row.frames + frames
                row.observed_dist = ((row.observed_dist*row.frames) + (observed_dist*frames))/new_frames
                row.observed_sig_pow = ((row.observed_sig_pow*row.frames) + (avg_sig_pow*frames))/new_frames
                row.observed_noise_pow = ((row.observed_noise_pow*row.frames) + (avg_noise_pow*frames))/new_frames
                row.observed_snr = value_in_db(row.observed_sig_pow/row.observed_noise_pow)
                row.frames = new_frames
                found = true
                break
            end
        end
    end
    if !found
        push!(df, [id, shap_rate, target_dist, observed_dist, target_snr, avg_sig_pow, avg_noise_pow, value_in_db(avg_sig_pow/avg_noise_pow), frames], promote=true)
    end
    CSV.write(channel_dm_path, df)
    return id
end

function channel_code_dump_bin_results(
    code::Code{BIT,BIT,REAL},
    constellation::AbstractString,
    snr::REAL,
    bec::NAT,
    fec::NAT,
    observed_dist::REAL,
    avg_sig_pow::REAL,
    avg_noise_pow::REAL,
    frames::NAT,
    change_chann::Bool
)
    @LOGD CHANNEL_SIM_DATA_TAG "dumping results..."
    constellation in ["OOK", "BPSK"] || ERROR(CHANNEL_SIM_DATA_TAG, "$constellation is not a known constellation")
    id = channel_code_dump_params(code, constellation,change_chann)
    channel_code_dump_sim_res(id, snr, bec, fec, frames, messagelength(code))
    channel_code_dump_bin_shap_res(id, shaping_rate(code), target_dist(code), snr, observed_dist, avg_sig_pow, avg_noise_pow, frames)
    @LOGD CHANNEL_SIM_DATA_TAG "dumping results done"
    return id
end