using CSV, DataFrames

const entropy_path_param = "src/channel_coding/polar_code/bit_channel_entropies/data/polarized_entropies_params.csv"
const entropy_path_data = "src/channel_coding/polar_code/bit_channel_entropies/data/polarized_entropies.csv"

function entropy_data_already_present(N::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    df = CSV.read(entropy_path_param, DataFrame)
    for row in eachrow(df)
        "$(row.blocklength)" == "$N" && "$(row.source_p)" == "$p1" && "$(row.channel)" == "$W" && return row.id
    end
    return -1
end

function entropy_fetch_hs(N::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    id = entropy_data_already_present(N, p1, W)
    id <= 0 && ERROR(FETCH_ENTROPY_TAG, "Not enough data for ($N, $p1, $W)")
    df = CSV.read(entropy_path_data, DataFrame)
    ready_in = zeros(Bool, N)
    vals = zeros(REAL, N)
    for row in eachrow(df[(df.id .== id),:])
        index = NAT(row.index)
        ready_in[index] && ERROR(FETCH_ENTROPY_TAG, "index $index present twicefor ($N, $p1, $W) ")
        vals[index] = row.value
        ready_in[index] = true
    end
    all(ready_in) || ERROR(FETCH_ENTROPY_TAG, "some index is missing in the data for ($N, $p1, $W)")
    return vals
end

function entropy_dump_params(N::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    id = entropy_data_already_present(N,p1,W)
    id > 0 && return id

    df = CSV.read(entropy_path_param, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, N, p1, "$W"], promote=true)
    return id, df
end

function entropy_id_already_used(id::NAT, df::DataFrame)
    any(df.id .== id)
end

function entropy_dump_results(id::NAT, hs::Vector{REAL}, frames::NAT)
    df = CSV.read(entropy_path_data, DataFrame)
    if entropy_id_already_used(id, df)
        used_index = zeros(Bool, length(hs))
        for row in eachrow(df[(df.id .== id),:])
            if row.id == id #just to be sure
                index = NAT(row.index)
                used_index[index] && ERROR(FETCH_ENTROPY_TAG, "Something wrong with incorportating data at id $id")
                row.value = (row.value*row.frames + hs[index]*frames) / (frames + row.frames)
                row.frames += frames
                used_index[index] = true
            end
        end
        all(used_index) || ERROR(FETCH_ENTROPY_TAG, "Something wrong with incorportating data at id $id")
    else
        for i in eachindex(hs)
            push!(df, [id, i, hs[i], frames], promote=true)
        end
    end
    return df
end

function entropy_dump_results(N::NAT, p1::REAL, W::Channel{BIT, Y}, hs::Vector{REAL}, frames::NAT) where {Y}
    id, df1 = entropy_dump_params(N, p1, W)
    df2 = entropy_dump_results(id, hs, frames)
    @LOG FETCH_ENTROPY_TAG "About to write CSV files 👀"
    CSV.write(entropy_path_data, df2)
    CSV.write(entropy_path_param, df1)
    @LOG FETCH_ENTROPY_TAG "Writing CSV files done 😮‍💨"
end

