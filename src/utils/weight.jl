function weight(i::NAT, n::NAT)
    acc = 0
    for j ∈ 1:n
        acc += i&1
        i >>>= 1
    end
    return acc
end
function weight(i::NAT)
    acc = 0
    while i != 0
        acc += i&1
        i >>>= 1
    end
    return acc
end
function weight(v::Vector{Bool})
    acc = 0
    for b ∈ v
        if b 
            acc += 1
        end
    end
    return acc
end
function weight(v::Vector{<:Integer})
    acc = 0
    for b ∈ v
        if b == 1
            acc += 1
        end
    end 
    return acc
end