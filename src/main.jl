include("../conf.jl")
include("utils/db.jl")
include("utils/weight.jl")
include("utils/entropy.jl")
include("utils/dic_parser.jl")
include("utils/opt_ook.jl")
include("utils/vec_gen.jl")
include("channel/channel.jl")
include("channel_coding/channel_sim.jl")
include("source_coding/lossy/lossy_coding_sim.jl")
include("distribution_matching/distribution_matching_sim.jl")

intereseting_params = Dict{Any, Any}(
    352 => 
        Dict{String, Any}(
            "rate" => 0.4,
            "decoding_list" => 32,
            "frozen_crc_width" => 5,
            "encoding_list" => 8,
            "info_k_pc_constraints" => 6,
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01],
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01],
            "shaping_k_pc_constraints" => 4,
            "shaping_rate" => 0.7319210675623926,
            "blocklength" => 128,
            "frozen_known_at_dec" => true,
            "dynamic_shaping_constraint" => "paritycheck",
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3,
            "shaping_known_at_dec" => true,
            "uniform_forzen_vec" => true,
            "dynamic_frozen_constraint" => "crc",
            "p" => 0.23962797516744064
        ), 
    312 => 
        Dict{String, Any}(
            "rate" => 0.4,
            "decoding_list" => 32,
            "frozen_crc_width" => 5,
            "encoding_list" => 8,
            "info_k_pc_constraints" => 6,
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01],
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01],
            "shaping_k_pc_constraints" => 6,
            "shaping_rate" => 0.7319210675623926,
            "blocklength" => 128,
            "frozen_known_at_dec" => true,
            "dynamic_shaping_constraint" => "rowmerge",
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    315 => 
        Dict{String, Any}(
            "rate" => 0.4,
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01],
            "shaping_k_pc_constraints" => 6,
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "rowmerge", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ), 
    374 => 
        Dict{String, Any}(
            "rate" => 0.4,
            "decoding_list" => 8,
            "frozen_crc_width" => 3, 
            "encoding_list" => 32,
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ),
    373 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 8, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ), 
    361 => 
        Dict{String, Any}(
            "rate" => 0.4,
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, "
            dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    368 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 8, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ), 
    380 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01],
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7319210675623926,
            "blocklength" => 128,
            "frozen_known_at_dec" => true,
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ),
    225 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "none", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    353 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 8, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    357 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 8, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    365 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7944210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        ),
    372 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 8, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ),
    309 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 1, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 6, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "rowmerge", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 3, 
            "shaping_known_at_dec" => true,
            "uniform_forzen_vec" => true,
            "dynamic_frozen_constraint" => "crc",
            "p" => 0.23962797516744064
        ),
    376 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "rowmerge", 
            "p" => 0.23962797516744064
        ),
    360 => 
        Dict{String, Any}(
            "rate" => 0.4, 
            "decoding_list" => 32, 
            "frozen_crc_width" => 5, 
            "encoding_list" => 32, 
            "info_k_pc_constraints" => 6, 
            "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
            "shaping_k_pc_constraints" => 4, 
            "shaping_rate" => 0.7319210675623926, 
            "blocklength" => 128, 
            "frozen_known_at_dec" => true, 
            "dynamic_shaping_constraint" => "paritycheck", 
            "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
            "shaping_crc_width" => 5, 
            "shaping_known_at_dec" => true, 
            "uniform_forzen_vec" => true, 
            "dynamic_frozen_constraint" => "crc", 
            "p" => 0.23962797516744064
        )
    )

intereseting_params =  Dict{Any, Any}(
        279 => 
            Dict{String, Any}(
                "rate" => 0.4,
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 6, 
                "shaping_rate" => 0.7944210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "crc", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "crc", 
                "p" => 0.23962797516744064
            ), 
       360 => 
            Dict{String, Any}(
                "rate" => 0.4, 
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 4, 
                "shaping_rate" => 0.7319210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "paritycheck", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "crc", 
                "p" => 0.23962797516744064
            ),
        #=
        353 => 
            Dict{String, Any}(
                "rate" => 0.4, 
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 8, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 4, 
                "shaping_rate" => 0.7944210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "paritycheck", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 3, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "crc", 
                "p" => 0.23962797516744064
            ),
        =#
        376 => 
            Dict{String, Any}(
                "rate" => 0.4, 
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 4, 
                "shaping_rate" => 0.7319210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "paritycheck", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "rowmerge", 
                "p" => 0.23962797516744064
            ),
        315 => 
            Dict{String, Any}(
                "rate" => 0.4,
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01],
                "shaping_k_pc_constraints" => 6,
                "shaping_rate" => 0.7319210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "rowmerge", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "crc", 
                "p" => 0.23962797516744064
            ), 
        225 => 
            Dict{String, Any}(
                "rate" => 0.4, 
                "decoding_list" => 32, 
                "frozen_crc_width" => 5, 
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 6, 
                "shaping_rate" => 0.7944210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "none", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "crc", 
                "p" => 0.23962797516744064
            ),
        #=
        296 => 
            Dict{String, Any}(
                "rate" => 0.4,
                "decoding_list" => 32,
                "frozen_crc_width" => 5,
                "encoding_list" => 32, 
                "info_k_pc_constraints" => 6, 
                "shaping_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "frozen_poly" => UInt8[0x01, 0x00, 0x01, 0x00, 0x00, 0x01], 
                "shaping_k_pc_constraints" => 6, 
                "shaping_rate" => 0.7319210675623926, 
                "blocklength" => 128, 
                "frozen_known_at_dec" => true, 
                "dynamic_shaping_constraint" => "crc", 
                "design_channel" => OOK(128, 0.7943282347242815, 0.2396279751674407),
                "shaping_crc_width" => 5, 
                "shaping_known_at_dec" => true, 
                "uniform_forzen_vec" => true, 
                "dynamic_frozen_constraint" => "rowmerge", 
                "p" => 0.23962797516744064
            ),
        =
        382 => 
            Dict{String,Any}(
                "blocklength" => 128,
                "design_channel" => OOK(128, sigma, optp),
                "p" => optp,
                "dynamic_shaping_constraint" => "none",
                "shaping_rate" => 0.6953125,
                "rate" => 0.4,
                "dynamic_frozen_constraint" => "crc",
                "frozen_crc_width" => 5,
                "frozen_poly" => get_default_poly(5),
                "frozen_known_at_dec" => true,
                "uniform_forzen_vec" => true,
                "encoding_list" => 32,
                "decoding_list" => 32,
            )
        =#
)

function general_sim0()
    while true 
        for dic in values(intereseting_params)
            code = get_polar_code(dic)
            id = channel_code_already_present(code, "OOK", true)
            if id ∉ keys(intereseting_params)
                @LOG MAIN_TAG "problem problem $id"
            end
            sim_bawgn(code, "OOK", [i/4 for i in 20:28], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true)
        end
    end
end

function general_sim()
    found = Dict()
    @LOG MAIN_TAG "HERE WE Go"
    for design_snr in [2.0]
        for blocklength in [128]
            for shap_const in ["none"]
                for froz_const in ["rowmerge"]
                    for encodingL in [32]
                        for info_k_pc_constraints in [4]
                            for shaping_k_pc_constraint in [4]
                                l = convert(NAT, log2(encodingL))
                                for decodingL in [32]
                                    L = convert(NAT, log2(decodingL))
                                    for unifo_fv in [true]
                                        sigma = snr_db_to_sigma(design_snr)
                                        optp, _ = opt_ook(sigma)
                                        h = h2(optp)
                                        shap_rates = [h + i/16 for i in -1:0]
                                        for shap_rate in shap_rates
                                            dic = Dict(
                                                "blocklength" => blocklength,
                                                "design_channel" => OOK(blocklength, sigma, optp),
                                                "p" => optp,
                                                "dynamic_shaping_constraint" => shap_const,
                                                "shaping_rate" => shap_rate,
                                                "rate" => 0.4,
                                                "shaping_crc_width" => max(l,3),
                                                "shaping_poly" => get_default_poly(max(l,3)),
                                                "shaping_known_at_dec" => true,
                                                "dynamic_frozen_constraint" => froz_const,
                                                "frozen_crc_width" => max(L,3),
                                                "frozen_poly" => get_default_poly(max(L,3)),
                                                "frozen_known_at_dec" => true,
                                                "uniform_forzen_vec" => unifo_fv,
                                                "encoding_list" => encodingL,
                                                "decoding_list" => decodingL,
                                                "info_k_pc_constraints" => info_k_pc_constraints,
                                                "shaping_k_pc_constraints" => shaping_k_pc_constraint
                                            )
                                            try
                                                code = get_polar_code(dic)
                                                id = channel_code_already_present(code, "OOK", true)
                                                if id == 242
                                                    println(dic)
                                                    break
                                                else
                                                    @LOG MAIN_TAG "Nope"
                                                end
                                                #sim_bawgn(code, "OOK", [i/4 for i in 16:32], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true)
                                            catch e
                                                
                                                @LOG MAIN_TAG "$e"
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    print(found)
end

function sim_row_merge()
    for decoding_list in [4, 8]
        dic = Dict(
            "blocklength" => 128,
            "design_channel" => BPSK(128, snr_db_to_sigma(2.0)),
            "p" => 0.5,
            "dynamic_shaping_constraint" => "none",
            "shaping_rate" => 0.0,
            "rate" => 60/128,
            "dynamic_frozen_constraint" => "rowmerge",
            "uniform_forzen_vec" => false,
            "encoding_list" => 1,
            "decoding_list" => decoding_list
        )
        code = get_polar_code(dic)
        sim_bawgn(code, "BPSK", [i/4 for i in 8:24], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true)

        dic = Dict(
            "blocklength" => 256,
            "design_channel" => BPSK(256, snr_db_to_sigma(2.0)),
            "p" => 0.5,
            "dynamic_shaping_constraint" => "none",
            "shaping_rate" => 0.0,
            "rate" => 75/256,
            "dynamic_frozen_constraint" => "rowmerge",
            "uniform_forzen_vec" => false,
            "encoding_list" => 1,
            "decoding_list" => decoding_list
        )
        code = get_polar_code(dic)
        sim_bawgn(code, "BPSK", [i/4 for i in 8:24], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true)
    end
end

function sim_vingart()
    optp, _ = opt_ook(snr_db_to_sigma(6.0))
    dic = Dict(
        "blocklength" => 1024,
        "design_channel" => OOK(1024, snr_db_to_sigma(6.0), 0.5),
        "p" => 0.5,
        "dynamic_shaping_constraint" => "none",
        "shaping_rate" => 0.0,
        "rate" => 0.67,
        "dynamic_frozen_constraint" => "crc",
        "frozen_crc_width" => 16,
        "frozen_poly" => get_default_poly(16),
        "frozen_known_at_dec" => true,
        "uniform_forzen_vec" => false,
        "encoding_list" => 1,
        "decoding_list" => 32
    )
    code = get_polar_code(dic)
    sim_bawgn(code, "OOK", [i/4 for i in 20:26], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS ÷ 2, DEFAULT_BATCH_SIZE, true)
    dic = Dict(
        "blocklength" => 1024,
        "design_channel" => OOK(1024, snr_db_to_sigma(6.0), optp),
        "p" => optp,
        "dynamic_shaping_constraint" => "none",
        "shaping_rate" => 115/1024,
        "rate" => 0.67,
        "dynamic_frozen_constraint" => "crc",
        "frozen_crc_width" => 16,
        "frozen_known_at_dec" => true,
        "frozen_poly" => get_default_poly(16),
        "uniform_forzen_vec" => true,
        "encoding_list" => 32,
        "decoding_list" => 32
    )
    code = get_polar_code(dic)
    sim_bawgn(code, "OOK", [i/4 for i in 20:26], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS ÷ 2, DEFAULT_BATCH_SIZE, true)
end

function sim_low_shaping()
    for change_chann in [true, false]
        for frozen_know in [true, false]
            #for d in reverse([0.0, 0.03125, 0.0625, 0.09375, 0.125, 0.15625, 0.1875, 0.21875, 0.25, 0.28125, 0.3125, 0.34375, 0.375, 0.40625])
            for d in reverse([0.0, 0.0625, 0.125, 0.25, 0.3125, 0.375, 0.40625, 0.125])    
                opt_p, _ = opt_ook(snr_db_to_sigma(2.0))
                dic = Dict(
                    "blocklength" => 32,
                    "design_channel" => OOK(32, snr_db_to_sigma(2.0), opt_p),
                    "p" => opt_p,
                    "dynamic_shaping_constraint" => "none",
                    "shaping_rate" => 1-d,
                    "rate" => 0.5,
                    "dynamic_frozen_constraint" => "crc",
                    "frozen_crc_width" => 3,
                    "frozen_known_at_dec" => frozen_know,
                    "frozen_poly" => BIT.([1, 1, 1, 0]),
                    "uniform_forzen_vec" => true,
                    "encoding_list" => 32,
                    "decoding_list" => 32
                )
                code = get_polar_code(dic)
                sim_bawgn(code, "OOK", [i/4 for i in 16:28], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS ÷ 2, DEFAULT_BATCH_SIZE, change_chann)
            end
        end
    end
end

function gen_5g_comp()
    for N in [64, 128, 256, 512]
        dic = Dict(
            "blocklength" => N,
            "design_channel" => BPSK(N, snr_db_to_sigma(DESIGN_SNR)),
            "p" => 0.5,
            "dynamic_shaping_constraint" => "none",
            "shaping_rate" => 1.0,
            "rate" => 0.5,
            "dynamic_frozen_constraint" => "crc",
            "frozen_crc_width" => 11,
            "frozen_known_at_dec" => false,
            "frozen_poly" => BIT.([1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
            "uniform_forzen_vec" => false,
            "encoding_list" => 1,
            "decoding_list" => 32
        )
        code = get_polar_code(dic)
        sim_bawgn(code, "BPSK", [i/4 for i in -4:16], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS ÷ 2, DEFAULT_BATCH_SIZE, true)
    end
end

function gen_5g_comp2()
    for N in [32]
        dic = Dict(
            "blocklength" => N,
            "design_channel" => BPSK(N, snr_db_to_sigma(DESIGN_SNR)),
            "p" => 0.5,
            "dynamic_shaping_constraint" => "none",
            "shaping_rate" => 1.0,
            "rate" => 0.5,
            "dynamic_frozen_constraint" => "crc",
            "frozen_crc_width" => 11,
            "frozen_known_at_dec" => false,
            "frozen_poly" => BIT.([1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
            "uniform_forzen_vec" => false,
            "encoding_list" => 1,
            "decoding_list" => 32
        )
        code = get_polar_code(dic)
        sim_bawgn(code, "BPSK", [i/2 for i in 14:14], 10, 10, 10, true)
    end
end

function gen_trad_lossy(p::REAL, encodingL::NAT, decodingL::NAT)
    for N in [64, 128, 256, 512]
        for targ_dist in [i/8 for i in 1:7]
            r = h2(p) - h2(targ_dist)
            dic = Dict(
                "blocklength" => N,
                "design_channel" => BSC(targ_dist),
                "p" => p,
                "shaping_rate" => h2(p),
                "rate" => r,
                "dynamic_frozen_constraint" => "none",
                "dynamic_shaping_constraint" => "none",
                "uniform_forzen_vec" => true,
                "encoding_list" => encodingL,
                "decoding_list" => decodingL
            )
            code = LossyPolarCode(get_polar_code(dic))
            lossy_coding_simulate_par(code, targ_dist, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE,true)
        end
    end
end

function general_sim2()

    for design_snr in [2.0]
        for blocklength in [64]
            for shap_const in ["paritycheck"]#["rowmerge", "crc", "paritycheck", "none"]
                for froz_const in ["none"]
                    for encodingL in [32]
                        for decodingL in [32]
                            for shaping_k_pc_constraint in [4,6,8,10,16,32,50]
                                for unifo_fv in [true]
                                    sigma = snr_db_to_sigma(design_snr)
                                    optp, _ = opt_ook(sigma)
                                    h = h2(optp)
                                    shap_rates = [h + i/16 for i in -10:1]
                                    for shap_rate in shap_rates
                                        dic = Dict(
                                            "blocklength" => blocklength,
                                            "design_channel" => OOK(blocklength, sigma, optp),
                                            "p" => optp,
                                            "dynamic_shaping_constraint" => shap_const,
                                            "shaping_rate" => shap_rate,
                                            "rate" => 0.5,
                                            "shaping_crc_width" => 5,
                                            "shaping_poly" => get_default_poly(5),
                                            "shaping_known_at_dec" => true,
                                            "dynamic_frozen_constraint" => froz_const,
                                            "frozen_crc_width" => 5,
                                            "frozen_poly" => get_default_poly(5),
                                            "frozen_known_at_dec" => true,
                                            "uniform_forzen_vec" => unifo_fv,
                                            "encoding_list" => encodingL,
                                            "decoding_list" => decodingL,
                                            "info_k_pc_constraints" => 6,
                                            "shaping_k_pc_constraints" => shaping_k_pc_constraint
                                        )
                                        try
                                            code = get_polar_code(dic)
                                            channel_code_already_present(code, "OOK", true) > 0 && continue
                                            sim_bawgn(code, "OOK", [i/2 for i in 4:12], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true)
                                        catch e
                                            println("main | $e")
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function try_shaping_with_channel()
    #conclusion : NOPE
    design_snr = 2.0
    sigma = snr_db_to_sigma(design_snr)
    optp, _ = opt_ook(sigma)
    for N in [128]
        dic = Dict(
            "blocklength" => N,
            "design_channel" => OOK(N, sigma, optp),
            "p" => optp,
            "uniform_forzen_vec" => true,
            "encoding_list" => 32,
            "decoding_list" => 32
        )
        code = create_good_ook_polar_code(dic)
        sim_bawgn(code, "OOK", [i/2 for i in 0:12], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true, false)
    end
end

function gen_none()
    design_snr = 2.0
    sigma = snr_db_to_sigma(design_snr)
    optp, _ = opt_ook(sigma)
    dic = Dict(
        "blocklength" => 128,
        "design_channel" => OOK(128, sigma, optp),
        "p" => optp,
        "dynamic_shaping_constraint" => "none",
        "shaping_rate" => 0.6953125,
        "rate" => 0.4,
        "dynamic_frozen_constraint" => "crc",
        "frozen_crc_width" => 5,
        "frozen_poly" => get_default_poly(5),
        "frozen_known_at_dec" => true,
        "uniform_forzen_vec" => true,
        "encoding_list" => 32,
        "decoding_list" => 32,
    )
    code = get_polar_code(dic)
    sim_bawgn(code, "OOK", [i/2 for i in 4:12], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true, true)
end

function sim_crc_prec()
    dic = intereseting_params[279]
    code = get_polar_code(dic)
    sim_bawgn(code, "OOK", [6.0], 1, 1, 1, true, false)
end

function lossy_coding_main()
    for blocklength in [128]
        for froz_const in ["rowmerge"]
            for shap_const in ["none"]
                for encodingL in [1, 8, 32]
                    for decodingL in [1, 8, 32]
                        encodingL != decodingL && continue
                        for info_k_pc_constraints in [4]
                            for shaping_k_pc_constraint in [4]
                                l = convert(NAT, log2(encodingL))
                                L = convert(NAT, log2(decodingL))
                                for unifo_fv in [false]
                                    for r in [i/16 for i in 1:15]
                                        dist = h2_inv(1-r)
                                        W = BSC(dist)
                                        dic = Dict(
                                            "blocklength" => blocklength,
                                            "design_channel" => W,
                                            "p" => 0.5,
                                            "dynamic_shaping_constraint" => shap_const,
                                            "shaping_rate" => 1.0,
                                            "rate" => r,
                                            "shaping_crc_width" => max(l,3),
                                            "shaping_poly" => get_default_poly(max(l,3)),
                                            "shaping_known_at_dec" => true,
                                            "dynamic_frozen_constraint" => froz_const,
                                            "frozen_crc_width" => max(L,3),
                                            "frozen_poly" => get_default_poly(max(L,3)),
                                            "frozen_known_at_dec" => true,
                                            "uniform_forzen_vec" => unifo_fv,
                                            "encoding_list" => encodingL,
                                            "decoding_list" => decodingL,
                                            "info_k_pc_constraints" => info_k_pc_constraints,
                                            "shaping_k_pc_constraints" => shaping_k_pc_constraint
                                        )
                                        try
                                            code = get_lossy_polar_code(dic)
                                            already_simulated(code) > 0 && continue
                                            lossy_coding_simulate_par(code, dist, 100000, DEFAULT_BATCH_SIZE)
                                        catch e
                                            @LOG MAIN_TAG "$e"
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function distribution_matching_main()
    while true 
        for blocklength in [128]
            for shap_const in ["paritycheck"]#["none", "crc", "paritycheck", "rowmerge"]
                for encodingL in [1]#[1, 8, 32]
                    for rr in [false]
                        if rr && encodingL != 1
                            @LOGD MAIN_TAG "Cannot use rr and encodingL != 1"
                            continue
                        else 
                            @LOGD MAIN_TAG "Let's do something else"
                        end
                        for shaping_k_pc_constraint in [4]
                            if shap_const != "paritycheck" && shaping_k_pc_constraint != 4 
                                @LOGD MAIN_TAG "Continue"
                                continue
                            end
                            p, _ = opt_ook(snr_db_to_sigma(2.0))
                            h = h2(p)
                            for r in [h + i/16 for i in -4:0]
                                @LOG MAIN_TAG "info rate is $r"
                                dic = Dict(
                                    "blocklength" => blocklength,
                                    "p" => p,
                                    "dynamic_shaping_constraint" => shap_const,
                                    "shaping_rate" => r,
                                    "shaping_crc_width" => 5,
                                    "shaping_poly" => get_default_poly(5),
                                    "shaping_known_at_dec" => true,
                                    "encoding_list" => encodingL,
                                    "shaping_k_pc_constraints" => shaping_k_pc_constraint,
                                    "randomized_rounding" => rr
                                )
                                try
                                    code = get_polar_code_distribution_matcher(dic)
                                    if already_simulated(code) > 0
                                        @LOG MAIN_TAG "already simulated, pass"
                                        continue
                                    end
                                    distribution_matching_simulate_par(code, 100000, DEFAULT_BATCH_SIZE)
                                catch e
                                    @LOG MAIN_TAG "$e"
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function distribution_matching_asymptot()
    for n in 11:15
        p, _ = opt_ook(snr_db_to_sigma(2.0))
        h = h2(p)
        r = h - 1/16
        N = 1 << n
        dic = Dict(
            "blocklength" => N,
            "p" => p,
            "dynamic_shaping_constraint" => "none",
            "shaping_rate" => r,
            "encoding_list" => 1,
            "randomized_rounding" => true
        )
        code = get_polar_code_distribution_matcher(dic)
        distribution_matching_simulate_par(code, 100000, DEFAULT_BATCH_SIZE)
    end
end

function crc_vs_row_merge()
    for N in [128]#, 256]
        @LOG MAIN_TAG "blocklength N=$N"
        sigma = snr_db_to_sigma(2.0)
        p, i = opt_ook(sigma)
        @LOG MAIN_TAG "found opt p $p and capacity $i"
        sr = h2(p)
        poly = N == 256 ? 
            BIT.([1,1,1,0,1,0,1,1,1,1,1]) :
            BIT.([1,1,1,0,0,1,0,1])
        for r in [0.3]
            @LOG MAIN_TAG "rate r $r"
            for sr in [h2(p) + j/16 for j in -1:-1]
                @LOG MAIN_TAG "shaping rate $sr"
                for froz_const in ["rowmerge", "crc"]
                    @LOG MAIN_TAG "frozen const $froz_const"
                    for encodingL in [32]
                        @LOG MAIN_TAG "encoding list $encodingL"
                        for decodinL in [32]#[8,32]
                            @LOG MAIN_TAG "getting dic"
                            dic = Dict(
                                "blocklength" => N,
                                "p" => p,
                                "dynamic_shaping_constraint" => "none",
                                "rate" => r,
                                "shaping_rate" => sr,
                                "frozen_crc_width" => N == 256 ? 10 : 7,
                                "frozen_poly" => poly,
                                "dynamic_frozen_constraint" => froz_const,
                                "design_channel" => OOK(N, sigma, p),
                                "frozen_known_at_dec" => true,
                                "uniform_forzen_vec" => true,
                                "encoding_list" => encodingL,
                                "decoding_list" => decodinL
                            )
                            try
                                @LOG MAIN_TAG "getting code"
                                code = get_polar_code(dic)
                                @LOG MAIN_TAG "simulating"
                                sim_bawgn(code, "OOK", [i/4 for i in 8:32], DEFAULT_MIN_ITERS, DEFAULT_MAX_ITERS, DEFAULT_BATCH_SIZE, true, true)
                            catch e
                                @LOG MAIN_TAG "e"
                            end
                        end
                    end
                end
            end
        end 
    end
end

function found_merges(N::NAT=128)
    sigma = snr_db_to_sigma(2.0)
    p, _ = opt_ook(sigma)
    W = OOK(N, sigma, p)
    hs, hw = get_entropies_for_source_and_channel(N, p, W)
    source_rel = sort(reverse(1:N), by=i->hs[i])
    print("r,sr,merges")
    for r in [i/32 for i in 1:31]
        for rs in [i/32 for i in 1:31]
            try
                D = source_rel[1:ceil(NAT, N*(1-rs))]
                I = target_next_rm_code(N, floor(NAT, N*r), hw, D)
                merges = merge_rows(I, N, D)
                println("$r,$rs,$(length(merges))")
            catch e
                @LOG MAIN_TAG "$e"
            end
        end
    end
end


#diag_weight_transform(64, 32)
#general_sim0()
#lossy_coding_main()
#distribution_matching_main()
#general_sim()
#distribution_matching_asymptot()
#general_sim()
#sim_crc_prec()
found_merges(128)
#crc_vs_row_merge()