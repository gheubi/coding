#=
    Convert a value from linear scale to decibels (dB) and vice versa.
    The functions value_in_db and value_from_db are used to convert a value v in linear scale to dB and back.
    The conversion is based on the formula: v_db = 10 * log10(v) for linear to dB, and v = 10^(v_db / 10) for dB to linear.
    @param v the value in linear scale
    @return the value in dB
=#
value_in_db(v::REAL) = 10*log10(v)
#=
    Convert a value from decibels (dB) to linear scale.
    The function value_from_db is used to convert a value v_db in dB to linear scale.
    The conversion is based on the formula: v = 10^(v_db / 10).
    @param v_db the value in dB
    @return the value in linear scale
=#
value_from_db(v_db::REAL) = 10^(v_db/10)