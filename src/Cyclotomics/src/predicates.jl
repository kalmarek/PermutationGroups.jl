function Base.hash(α::Cyclotomic, h::UInt)
    normalform!(α)
    return hash(coeffs(α), hash(conductor(α), hash(Cyclotomic, h)))
end

function Base.:(==)(α::Cyclotomic, β::Cyclotomic)
    if conductor(α) == conductor(β)
        normalform!(α)
        normalform!(β)
        return coeffs(α) == coeffs(β)
    else
        l = lcm(conductor(α), conductor(β))
        return embed(α, l) == embed(β, l)
    end
end

Base.iszero(α::Cyclotomic) =
    all(iszero, values(α)) || (normalform!(α); all(iszero, values(α)))

Base.isreal(α::Cyclotomic) = conductor(reduced_embedding(α)) == 1

function Base.isone(α::Cyclotomic)
    β = reduced_embedding(α)
    conductor(β) == 1 || return false
    return isone(β[0])
end

function isnormalized(
    α::Cyclotomic,
    basis = first(zumbroich_viacomplement(conductor(α))),
)
    # return all(in(basis), exponents(α))
    for e in exponents(α)
        e in basis || return false
    end
    return true
end

function _all_equal(α::Cyclotomic, exps, value = α[first(exps)])
    # return all(e -> α[e] == val, exps)
    for e in exps
        α[e] == value || return false
    end
    return true
end
