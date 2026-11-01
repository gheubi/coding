abstract type Channel{X,Y} end
transmit!(W::Channel{X,Y}, src::AbstractVector{X}, dst::AbstractVector{Y}) where {X, Y} = ABSTRACT(CHANNEL__TAG, W, "transmit!")
transmit_one_symbol(W::Channel{X,Y}, x::X) where {X,Y} = ABSTRACT(CHANNEL__TAG, "transmit_one_symbol")
(W::Channel{X,Y})(y::Y, x::X) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "transition probability")
LLR(W::Channel{BIT,Y}, y::Y) where {Y} = log(W(y,BIT(0))) - log(W(y,BIT(1)))

last_gen_noise_pow(W::Channel{X,Y}) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "last_gen_noise_pow")
last_signal_pow(W::Channel{X,Y}) where {X,Y} = ABSTRACT(CHANNEL__TAG, W, "last_signal_pow")

Base.show(io::IO, W::Channel) = ABSTRACT(CHANNEL__TAG, W, "Base.show")
Base.parse(::Type{T}, s::AbstractString) where {T<:Channel} = ABSTRACT(CHANNEL__TAG, T, "Base.parse()")

include("awgnc.jl")
include("binary_constellation.jl")
include("bec.jl")
include("bsc.jl")
include("noisless_channel.jl")