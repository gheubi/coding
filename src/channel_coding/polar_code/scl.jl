mutable struct PcWrapper
    pc::PolarCode
end

struct SCLRes
    paths::Vector{Vector{BIT}}    #of size (L, N)
    llr::Vector{Vector{REAL}}   #of size (L, N)
    metrics::Vector{REAL}       #of size (L,)
end

Base.show(io::IO, s::SCLRes) = begin
    st = ""
    for i in eachindex(s.paths)
        st = st * "[$i] ($(s.metrics[i])) | $(print_vector(s.paths[i]))\n"
    end
    print(io, st)
end

function emptySCL(N::NAT, L::NAT)
    paths = Vector{Vector{BIT}}(undef, L)
    llrs = Vector{Vector{REAL}}(undef, L)
    @inbounds for i in 1:L
        paths[i] = zeros(BIT, N)
        llrs[i] = zeros(REAL, N) 
    end
    SCLRes(paths, llrs, zeros(REAL, L))
end

mutable struct SCLResC
    paths::Ptr{Ptr{UInt8}}        # bit **paths
    llrs::Ptr{Ptr{Float64}}      # real **llrs (real = double)
    metrics::Ptr{Float64}        # real *metrics
    N::UInt16                    # nat N
    paths_size::UInt16           # nat paths_size
end
mutable struct SCLWorkingDS
    ptr::Ptr{Cvoid}
end
mutable struct SCLResDS
    ptr::Ptr{Cvoid}
end

function get_scl_working_ds(N::NAT, L::NAT)
    res = ccall(
        (:init_llrtreelist_working_args, "src/channel_coding/polar_code/llr_tree.so"),
        Ptr{Cvoid},
        (UInt16, UInt16),
        UInt16(N),
        UInt16(L)
    )::Ptr{Cvoid}
    working_args = SCLWorkingDS(res)
    finalizer(working_args) do wa
        ccall(
        (:deinit_llrtreelist_working_args, "src/channel_coding/polar_code/llr_tree.so"),
        Cvoid,
        (Ptr{Cvoid},),
        wa.ptr
        )
    end
    return working_args
end
function get_scl_res_ds(N::NAT, L::NAT)
    res = ccall(
        (:tree_list_init, "src/channel_coding/polar_code/llr_tree.so"),
        Ptr{Cvoid},
        (UInt16, UInt16),
        UInt16(N),
        UInt16(L)
    )::Ptr{Cvoid}
    res_args = SCLResDS(res)
    finalizer(res_args) do ra 
        ccall(
        (:tree_list_deinit, "src/channel_coding/polar_code/llr_tree.so"),
        Cvoid,
        (Ptr{Cvoid},),
        ra.ptr
        )
    end
    return res_args
end
function convert_to_scl_res!(c_result::SCLResC, dst::SCLRes)
    N = Int(c_result.N)
    N == length(dst.llr[1]) || 
        ERROR(PC_SCL_TAG, "resulting path length doesn't match expected ($(length(dst.llr[1])) != $N)")
    paths_size = Int(c_result.paths_size)
    paths_size == length(dst.metrics) || 
        ERROR(PC_SCL_TAG, "resulting number of path doesn't match l ($(paths_size) != $(length(dst.metrics)))")
    
    path_ptrs = unsafe_wrap(Vector{Ptr{UInt8}}, c_result.paths, paths_size)
    paths = dst.paths
    for i in 1:paths_size
        copy!(paths[i], unsafe_wrap(Vector{UInt8}, path_ptrs[i], N))
    end
    
    llrs_ptrs = unsafe_wrap(Vector{Ptr{Float64}}, c_result.llrs, paths_size)
    llrs = dst.llr
    for i in 1:paths_size
        copy!(llrs[i], unsafe_wrap(Vector{Float64}, llrs_ptrs[i], N))
    end
    
    copy!(dst.metrics, unsafe_wrap(Vector{Float64}, c_result.metrics, paths_size))
    
    return dst
end
function convert_to_scl_res!(res::SCLResDS, dst::SCLRes)
    c_result = unsafe_load(Ptr{SCLResC}(res.ptr))
    return convert_to_scl_res!(c_result, dst)
end
function c_decide_dynamical_frozen_bit(pc_ptr::Ptr{Cvoid}, u_ptr::Ptr{UInt8}, index::UInt16, N::UInt16)
    pc_ref = unsafe_pointer_to_objref(pc_ptr)::PcWrapper
    u = unsafe_wrap(Vector{UInt8}, u_ptr, N; own=false)
    index_conv = convert(NAT, index)
    return decide_dynmically_frozen_bit(pc_ref.pc, u, index_conv)
end
const c_dynamic_frozen_bit_at = @cfunction(
    c_decide_dynamical_frozen_bit,
    UInt8,
    (Ptr{Cvoid}, Ptr{UInt8}, UInt16, UInt16)
)

function scl_decode!(
    llrs_vec::Vector{REAL},
    L::NAT,
    pc::PolarCode,
    ptr_w::SCLWorkingDS,
    ptr_r::SCLResDS,
    statically_frozen_set::Vector{Bool},
    statically_frozen_vec::Vector{BIT},
    dynamically_frozen_set::Vector{Bool},
    dst::SCLRes
    )
    pc_wrap = PcWrapper(pc)
    result = GC.@preserve pc_wrap ptr_w ptr_r llrs_vec statically_frozen_set statically_frozen_vec dynamically_frozen_set begin
        pc_ptr = pointer_from_objref(pc_wrap)
        ccall(
            (:llrtreelist, "src/channel_coding/polar_code/llr_tree.so"),
            Int64,
            (Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Bool}, Ptr{UInt8}),
            llrs_vec,
            ptr_w.ptr,
            ptr_r.ptr,
            pc_ptr,
            statically_frozen_set,
            statically_frozen_vec,
            dynamically_frozen_set,
            c_dynamic_frozen_bit_at
        )
    end
    result >= 0 || ERROR(PC_SCL_TAG, "decoding failed in the C-code")
    return convert_to_scl_res!(ptr_r, dst)
end

function mostlikely!(scl::SCLRes, dst::Vector{BIT})
    isempty(scl.metrics) && ERROR(PC_SCL_TAG, "empty scl metrics")
    minI = argmin(i -> iszero(scl.metrics[i]) ? ∞ : scl.metrics[i], eachindex(scl.metrics))
    return copy!(dst, scl.paths[minI])
end

function mostlikely!(scl::SCLRes, pc::PolarCode, dst::Vector{BIT})
    metrics = scl.metrics
    isempty(metrics) && ERROR(PC_SCL_TAG, "empty scl metrics")
    first_index = 1
    best_metric = ∞
    best_index = 0
    @inbounds for i in eachindex(metrics)
        @LOGV PC_SCL_TAG "considering metric $i ($(metrics[i]))"
        m = metrics[i]
        if m == 0
            @LOGV PC_SCL_TAG "metric is 0, pass"
            continue
        end
        if m < metrics[first_index]
            @LOGV PC_SCL_TAG "metric improves best metric seen so far"
            first_index = i
        end
        isnan(m) && continue
        if m < best_metric && is_valid_message(pc, scl.paths[i])
            @LOGV PC_SCL_TAG "metric is of a valid path"

            best_metric = m
            best_index = i
        end
    end
    ptr = best_index != 0 ? scl.paths[best_index] : scl.paths[first_index]
    return copy!(dst, ptr)
end 