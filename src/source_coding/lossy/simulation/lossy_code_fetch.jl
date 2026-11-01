using CSV, DataFrames

const lossy_param_path = "src/source_coding/lossy/simulation/data/lossy_coding_params.csv"
const lossy_data_path = "src/source_coding/lossy/simulation/data/lossy_coding_sim.csv"

function already_simulated(code::LossyCode{U,X,Y}) where {U,X,Y}
    df = CSV.read(lossy_param_path, DataFrame)
    c = "$code"
    d = "$(input_distribution(code))"
    r = rate(code)
    for row in eachrow(df)
        @LOGV LOSSY_SIM_DATA_TAG "Already simulated code with id $(row.id)"
        row.code == c && row.source == d && row.rate == r && return row.id
    end
    return -1
end 

function dump_params(code::LossyCode{U,X,Y}) where {U,X,Y}
    id = already_simulated(code)
    id > 0 && return id
    df = CSV.read(lossy_param_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, rate(code), "$(input_distribution(code))", "$code"], promote=true)
    CSV.write(lossy_param_path, df)
    return id
end

function dump_lossy_res(id::NAT, rate::REAL, target_dist::REAL, dist::REAL, frames::NAT)
    df = CSV.read(lossy_data_path, DataFrame)
    found = false
    if any(df.id .== id)
        for row in eachrow(df)
            if row.id == row.id && row.rate == rate && row.target_distortion == target_dist
                row.dist = (row.dist*row.frames + dist*frames)/(row.frames + frames)
                row.frames += frames
                found = true
                break
            end
        end
    end
    if !found
        push!(df, [id, rate, target_dist, dist, frames], promote=true)
    end
    CSV.write(lossy_data_path, df)
    return id
end

function dump_lossy_sim_res(code::LossyCode{U,X,Y}, target_dist::REAL, dist::REAL, frames::NAT) where {U,X,Y}
    id = dump_params(code)
    dump_lossy_res(id, rate(code), target_dist, dist, frames)
    return id
end


