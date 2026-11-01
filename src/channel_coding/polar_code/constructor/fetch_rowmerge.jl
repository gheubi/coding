using CSV, DataFrames


const rm_params_path = "data/code/polar/row_merges_param.csv"
const rm_data_path = "data/code/polar/row_merges.csv"

function merges_already_present(blocklength::NAT, rate_profile::Vector{NAT}, forbiden_set::Vector{NAT}) 
    rp = "$rate_profile"
    fs = "$forbiden_set"
    df = CSV.read(rm_params_path, DataFrame)
    for row in eachrow(df)
        "$(row.blocklength)" == "$blocklength" &&  row.rate_profile == rp && row.forbiden_set == fs && return row.id
    end
    return -1
end

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

function merges_dump(id::NAT, merges::Vector{NTuple{2,NAT}})
    df = CSV.read(rm_data_path, DataFrame)
    for t in merges
        push!(df, [id, t[1], t[2]], promote=true)
    end
    CSV.write(rm_data_path, df)
    return id
end

function dump_merges_param(blocklength::NAT, rate_profile::Vector{NAT}, forbiden_set::Vector{NAT})
    id = merges_already_present(blocklength, rate_profile, forbiden_set)
    id > 0 && return id
    df = CSV.read(rm_params_path, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, blocklength, "$rate_profile", "$forbiden_set"], promote=true)
    CSV.write(rm_params_path, df)
    return id
end

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