using CSV, DataFrames


const rm_params_path = "data/code/polar/row_merges_param.csv"
const rm_data_path = "data/code/polar/row_merges.csv"

#=
    Check if the row merges for a given blocklength, rate profile and forbiden set are already present in the CSV file.
    The function reads the CSV file containing the parameters of previously computed row merges and checks if a set of merges with the same parameters exists.
    If such a set exists, it returns its ID; otherwise, it returns -1.
    @param blocklength the blocklength of the code
    @param rate_profile the rate profile of the code
    @param forbiden_set the forbiden set of the code
    @return the ID of the existing set of merges or -1 if not found
=#
function merges_already_present(blocklength::NAT, rate_profile::Vector{NAT}, forbiden_set::Vector{NAT}) 
    rp = "$rate_profile"
    fs = "$forbiden_set"
    df = CSV.read(rm_params_path, DataFrame)
    for row in eachrow(df)
        "$(row.blocklength)" == "$blocklength" &&  row.rate_profile == rp && row.forbiden_set == fs && return row.id
    end
    return -1
end

#=
    Fetch the row merges for a given ID from the CSV file.
    The function reads the CSV file containing the row merges and constructs a vector of merges corresponding to the given ID.
    @param id the ID of the set of merges
    @return a vector of row merges
=#
function merges_fetch(id::NAT)
    df = CSV.read(rm_data_path, DataFrame)
    merges = Vector{NTuple{2,NAT}}([])
    for row in eachrow(df)
        if row.id == id
            push!(merges, (row.i, row.f))
        end
    end
    return merges
end

#=
    Dump the row merges for a given ID to the CSV file.
    The function reads the CSV file containing the row merges and adds the new merges corresponding to the given ID.
    @param id the ID of the set of merges
    @param merges a vector of row merges
    @return the ID of the set of merges
=#
function merges_dump(id::NAT, merges::Vector{NTuple{2,NAT}})
    df = CSV.read(rm_data_path, DataFrame)
    for t in merges
        push!(df, [id, t[1], t[2]], promote=true)
    end
    CSV.write(rm_data_path, df)
    return id
end

#=
    Dump the parameters of a set of row merges to the CSV file.
    The function checks if the merges have already been computed and returns their ID if so.
    Otherwise, it adds the parameters to the CSV file and returns the new ID.
    @param blocklength the blocklength of the code
    @param rate_profile the rate profile of the code
    @param forbiden_set the forbiden set of the code
    @return the ID of the set of merges
=#

function dump_merges_param(blocklength::NAT, rate_profile::Vector{NAT}, forbiden_set::Vector{NAT})
    id = merges_already_present(blocklength, rate_profile, forbiden_set)
    id > 0 && return id
    df = CSV.read(rm_params_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, blocklength, "$rate_profile", "$forbiden_set"], promote=true)
    CSV.write(rm_params_path, df)
    return id
end

#=
    Dump the parameters and results of a set of row merges to the CSV files.
    The function first checks if the merges are already present using merges_already_present. If they are not present, 
    it dumps the parameters using dump_merges_param and then dumps the results using merges_dump.
    Finally, it returns the ID of the set of merges.
    @param blocklength the blocklength of the code
    @param rate_profile the rate profile of the code
    @param forbiden_set the forbiden set of the code
    @param merges a vector of row merges
    @return the ID of the set of merges
=#
function dump_merges(blocklength::NAT, rate_profile::Vector{NAT}, forbiden_set::Vector{NAT}, merges::Vector{NTuple{2,NAT}})
    id = merges_already_present(blocklength, rate_profile, forbiden_set)
    if id > 0
        merges2 = merges_fetch(id)
        for t in merges2
            t in merges || ERROR(ROWMERGE_DATA_TAG, "merges already in but don't match $merges2 != $merges")
        end
    else
        id = dump_merges_param(blocklength, rate_profile, forbiden_set)
        merges_dump(id, merges)
    end
    return id
end