using CSV, DataFrames

const lossy_param_path = "data/code/simulation/lossy_coding_params.csv"
const lossy_data_path = "data/code/simulation/lossy_coding_sim.csv"

#=
    Check if a lossy code with the given parameters has already been simulated.
    The function reads the CSV file containing the parameters of previously simulated lossy codes and checks if a code with the same parameters exists.
    If such a code exists, it returns its ID; otherwise, it returns -1.
    @param code the lossy code to check
    @return the ID of the existing code or -1 if not found
=#
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

#=
    Dump the parameters of a lossy code to the CSV file.
    The function checks if the code has already been simulated and returns its ID if so.
    Otherwise, it adds the code's parameters to the CSV file and returns the new ID.
    @param code the lossy code whose parameters to dump
    @return the ID of the code
=#
function dump_params(code::LossyCode{U,X,Y}) where {U,X,Y}
    id = already_simulated(code)
    id > 0 && return id
    df = CSV.read(lossy_param_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, rate(code), "$(input_distribution(code))", "$code"], promote=true)
    CSV.write(lossy_param_path, df)
    return id
end

#=
    Dump the results of a lossy code simulation to the CSV file.
    The function checks if the code has already been simulated and updates the results if so.
    Otherwise, it adds the new results to the CSV file.
    @param id the ID of the lossy code
    @param rate the rate of the lossy code
    @param target_dist the target distortion of the lossy code
    @param dist the achieved distortion of the lossy code
    @param frames the number of frames simulated
    @return the ID of the code
=#
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

#=
    Dump the results of a lossy code simulation to the CSV file.
    The function first dumps the parameters of the lossy code and then dumps the results.
    @param code the lossy code whose results to dump
    @param target_dist the target distortion of the lossy code
    @param dist the achieved distortion of the lossy code
    @param frames the number of frames simulated
    @return the ID of the code
=#
function dump_lossy_sim_res(code::LossyCode{U,X,Y}, target_dist::REAL, dist::REAL, frames::NAT) where {U,X,Y}
    id = dump_params(code)
    dump_lossy_res(id, rate(code), target_dist, dist, frames)
    return id
end


