const NAT = Int64
const BIT = UInt8
const REAL = Float64
const ∞ = Inf64
const ERASURE = 0x3
const CHUNK_T = UInt8
const CHUNK_L = 8
const DESIGN_SNR = 2.0
const NINTERACTIVE = Threads.nthreads(:interactive)
const NDEFAULT = Threads.nthreads(:default)
const LOG_ACT = false
const LOG_DEBUG_ACT = false
const LOG_DEBUG_ALLOWED_TAG = ["PolarCode", "Channel Simulation"]
const LOG_VERBOSE_ACT = false
const LOG_VERBOSE_ALLOWED_TAG = ["Successive Cancellation Decoder"]
const CHECK_ARGS = false 
const DEFAULT_MIN_ITERS = 100_000
const DEFAULT_MAX_ITERS = 1_000_000
const DEFAULT_BATCH_SIZE = 20_000

const MAIN_TAG = "Main"
const AWGNC_TAG = "AWGNC"
const BEC_TAG = "BEC"
const BIN_CONST_AWGNC_TAG = "Binary Constellation over AWGNC"
const BSC_TAG = "BSC"
const CHANNEL_TAG = "Channel"
const NOISELESS_CHANEL_TAG = "Noiseless Channel"
const FETCH_ENTROPY_TAG = "Fetch Polarized Entropies"
const POLARIZE_ENTROPY_TAG = "Polarizing Entropies"
const ROWMERGE_DATA_TAG = "RowMerge Data"
const POLAR_CODE_CONSTRUCTOR_TAG = "PolarCode Constructor"
const ROWMERGE_TAG = "RowMerge"
const PCDC_CONSTRAINT_CRC_TAG = "CRC Dynamic Constraint"
const PCDC_CONSTRAINT_TAG = "Dynamic Constraint"
const PCDC_CONSTRAINT_NONE_TAG = "None Dynamic Constraint"
const PCDC_CONSTRAINT_PC_TAG = "Parity Check Dynmaic Constraint"
const PCDC_CONSTRAINT_ROWMERGE_TAG = "RowMerge Dynamic Constraint"
const PC_UTILS_WEIGHT_ENUMERATOR = "PolarCode Weight Enumator"
const POLAR_CODE_TAG = "PolarCode" 
const POLAR_TRANSFORM = "Polar Transform"
const PC_SCL_TAG = "Successive Cancellation List Decoder"
const CHANNEL_SIM_DATA_TAG = "Channel Simulation Data"
const CHANNEL_CODE_TAG = "Channel Code"
const CHANNEL_SIM_TAG = "Channel Simulation"
const LOSSY_POLAR_CODE_TAG = "Lossy Polar Code"
const LOSSY_SIM_DATA_TAG = "Lossy Simulation Data"
const LOSSY_CODE_TAG = "Lossy Code"
const LOSSY_SIM_TAG = "Lossy Simulation"
const UTILS_COMPRESS_TAG = "Compressor Util"
const UTILS_CRC_TAG = "CRC Util : I'm deprecated, please don't use me 😢 => you can adapt the crc algo in $PCDC_CONSTRAINT_CRC_TAG"
const UTILS_DB_TAG = "dB Util"
const UTILS_PARSER_TAG = "Parser Util"
const UTILS_ENTROPY_TAG = "Entropy Util"
const UTILS_OOK_OPTIM_TAG = "On-Off Keying Optimzer Util"
const UTILS_VECT_GENERATOR_TAG = "Vector Generator Util"
const UTILS_WEIGHT = "Weight Util"
const DM_TAG = "Distribution Matcher"
const DM_SIM_TAG = "Distribution Matching Simulation"
const PCDM_TAG = "Polar Code Distribution Matcher"
const RR_SC_TAG = "Randomized Rounding Successive Cancellation"


function default_tid()
    if (Threads.threadpool() == :default)
        return Threads.threadid() - NINTERACTIVE
    else 
        return 1
    end
end

ERROR(TAG, msg) = error("$TAG | $msg 😒")
ABSTRACT(TAG, c, f_name) = ERROR(TAG, "$f_name not implemented for $(typeof(c)) 😒")
ABSTRACT(TAG, ::Type{T}, f_name) where {T} = ERROR(TAG, "$f_name not implemented for $T 😒")

macro LOG(TAG, msg)
    if LOG_ACT
        return quote
            println(string($TAG, " | ", $msg, " [", default_tid(), "]"))
        end |> esc
    else
        return :(nothing)
    end
end 

macro LOGD(TAG, msg)
    if LOG_DEBUG_ACT
        return quote
            println(string($TAG, " | ", $msg, " [", default_tid(), "]"))
        end |> esc
    else
        return :(nothing)
    end
end 

macro LOGV(TAG, msg)
    if LOG_VERBOSE_ACT
        return quote
            println(string($TAG, " | ", $msg, " [", default_tid(), "]"))
        end |> esc
    else
        return :(nothing)
    end
end 

macro check(ex)
    if CHECK_ARGS
        return esc(ex)
    else
        return nothing
    end
end
