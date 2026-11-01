using CSV, DataFrames

const entropy_path_param = "data/code/polar/polarized_entropies_params.csv" 
const entropy_path_data = "data/code/polar/polarized_entropies.csv"

#=
    Check if the polarized entropies for a given blocklength, source probability, and channel are already present in the CSV file.
    The function reads the CSV file containing the parameters of previously computed polarized entropies and checks if a set of entropies with the same parameters exists.
    If such a set exists, it returns its ID; otherwise, it returns -1.
    @param N the blocklength of the code
    @param p1 the source probability
    @param W the channel
    @return the ID of the existing set of entropies or -1 if not found
=#
function entropy_data_already_present(N::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    df = CSV.read(entropy_path_param, DataFrame)
    for row in eachrow(df)
        "$(row.blocklength)" == "$N" && "$(row.source_p)" == "$p1" && "$(row.channel)" == "$W" && return row.id
    end
    return -1
end

#=
    Fetch the polarized entropies for a given blocklength, source probability, and channel from the CSV file.
    The function first checks if the entropies are already present using entropy_data_already_present. If they are not present, it raises an error.
    If they are present, it reads the CSV file containing the entropies and constructs a vector of entropies corresponding to the given parameters.
    @param N the blocklength of the code
    @param p1 the source probability
    @param W the channel
    @return a vector of polarized entropies
=#
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

#=
    Dump the parameters of a set of polarized entropies to the CSV file.
    The function checks if the entropies have already been computed and returns their ID if so.
    Otherwise, it adds the parameters to the CSV file and returns the new ID.
    @param N the blocklength of the code
    @param p1 the source probability
    @param W the channel
    @return the ID of the set of entropies
=#
function entropy_dump_params(N::NAT, p1::REAL, W::Channel{BIT,Y}) where {Y}
    id = entropy_data_already_present(N,p1,W)
    id > 0 && return id

    df = CSV.read(entropy_path_param, DataFrame)
    id = length(df.id) == 0 ? 1 : maximum(df.id) + 1
    push!(df, [id, N, p1, "$W"], promote=true)
    return id, df
end

# cehck if the id is already used in the data file
function entropy_id_already_used(id::NAT, df::DataFrame)
    any(df.id .== id)
end

#=
    Dump the results of a set of polarized entropies to the CSV file.
    The function checks if the entropies have already been computed and updates the results if so.
    Otherwise, it adds the new results to the CSV file.
    @param id the ID of the set of entropies
    @param hs a vector of polarized entropies
    @param frames the number of frames simulated
    @return the updated DataFrame containing the dumped results
=#
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

#=
    Dump the parameters and results of a set of polarized entropies to the CSV files.
    The function first dumps the parameters using entropy_dump_params and then dumps the results using entropy_dump_results.
    Finally, it writes the updated DataFrames to their respective CSV files.
    @param N the blocklength of the code
    @param p1 the source probability
    @param W the channel
    @param hs a vector of polarized entropies
    @param frames the number of frames simulated
=#
function entropy_dump_results(N::NAT, p1::REAL, W::Channel{BIT, Y}, hs::Vector{REAL}, frames::NAT) where {Y}
    id, df1 = entropy_dump_params(N, p1, W)
    df2 = entropy_dump_results(id, hs, frames)
    @LOG FETCH_ENTROPY_TAG "About to write CSV files 👀"
    CSV.write(entropy_path_data, df2)
    CSV.write(entropy_path_param, df1)
    @LOG FETCH_ENTROPY_TAG "Writing CSV files done 😮‍💨"
end

