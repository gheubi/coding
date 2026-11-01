value_in_db(v::REAL) = 10*log10(v)
value_from_db(v_db::REAL) = 10^(v_db/10)