using CSV, DataFrames

const distribution_matching_params_path = "src/distribution_matching/simulation/data/distribution_matching_params.csv"
const distribution_matching_data_path = "src/distribution_matching/simulation/data/distribution_matching_sim.csv"

function already_simulated(dm::DistributionMatcher{X}) where {X}
    df = CSV.read(distribution_matching_params_path, DataFrame)
    d = "$dm"
    r = info_rate(dm)
    td = "$(target_distribution(dm))"
    bl = blocklength(dm)
    for row in eachrow(df)
        row.distribution_matcher == d && row.info_rate == r && 
            row.target_distribution == td && row.blocklength == bl && return row.id
    end
    return -1
end

function dump_params(dm::DistributionMatcher{X}) where {X}
    id = already_simulated(dm)
    id > 0 && return id
    df = CSV.read(distribution_matching_params_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, info_rate(dm), "$(target_distribution(dm))", blocklength(dm), "$dm"], promote=true)
    CSV.write(distribution_matching_params_path, df)
    return id
end

function dump_dm_res(id::NAT, average_distribution::REAL, fer::REAL, frames::NAT)
    df = CSV.read(distribution_matching_data_path, DataFrame)
    found = false
    if any(df.id .== id)
        for row in eachrow(df)
            if row.id == row.id
                row.average_distribution = (row.average_distribution*row.frames + average_distribution*frames) / (row.frames + frames)
                row.fer = (row.fer*row.frames + fer*frames) / (row.frames + frames)
                row.frames += frames
                found = true
                break
            end
        end
    end
    if !found
        push!(df, [id, average_distribution, fer, frames], promote=true)
    end
    CSV.write(distribution_matching_data_path, df)
    return id
end

function dump_distribution_matching_sim_res(dm::DistributionMatcher{X}, average_distribution::REAL, fer::REAL, frames::NAT) where {X}
    id = dump_params(dm)
    dump_dm_res(id, average_distribution, fer, frames)
    return id
end
