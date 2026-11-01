#=
    Abstract channel type and interface for all channels.
    The generic types X and Y are the input and output symbol types, respectively.
=#
abstract type Channel{X,Y} end

#=
    Transmit a vector of symbols through the channel.
    The input symbols are in `src` and the output symbols are stored in `dst`.
=#
transmit!(W::Channel{X,Y}, src::AbstractVector{X}, dst::AbstractVector{Y}) where {X, Y} = ABSTRACT(CHANNEL__TAG, W, "transmit!")

#=
    Transmit a single symbol through the channel.
    The input symbol is `x` and the output symbol is returned.
=#
transmit_one_symbol(W::Channel{X,Y}, x::X) where {X,Y} = ABSTRACT(CHANNEL__TAG, "transmit_one_symbol")

#=
    Get the transition probability of the channel.
    The input symbol is `x` and the output symbol is `y`.
=#
(W::Channel{X,Y})(y::Y, x::X) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "transition probability")

#=
    For binary input channels, get the log-likelihood ratio (LLR)
    of the channel output.
    The input symbol is `y` and the output is the LLR value.
=#
LLR(W::Channel{BIT,Y}, y::Y) where {Y} = log(W(y,BIT(0))) - log(W(y,BIT(1)))

#=
    Get the last generated noise power of the channel, which is the number of symbols that were flipped during the last 
    call to transmit!().
=#
last_gen_noise_pow(W::Channel{X,Y}) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "last_gen_noise_pow")
#=
    Get the last signal power of the channel, which is the number of symbols that were transmitted during the last 
    call to transmit!().
=#
last_signal_pow(W::Channel{X,Y}) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "last_signal_pow")

#=
    Show the channel information.
=#
Base.show(io::IO, W::Channel) = ABSTRACT(CHANNEL__TAG, W, "Base.show")
#=
    Parse a channel from a string.
=#
Base.parse(::Type{T}, s::AbstractString) where {T<:Channel} = ABSTRACT(CHANNEL__TAG, T, "Base.parse()")

include("awgnc.jl")
include("binary_constellation.jl")
include("bec.jl")
include("bsc.jl")
include("noisless_channel.jl")