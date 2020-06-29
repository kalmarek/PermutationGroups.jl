function Base.hash(α::Cyclotomic, h::UInt)
    normalform!(α)
    return hash(coeffs(α), hash(conductor(α), hash(Cyclotomic, h)))
end

Base.:(==)(α::Cyclotomic, β::Cyclotomic) =
    ==(α, β, Val{conductor(α) == conductor(β)}())

function Base.:(==)(α::Cyclotomic, β::Cyclotomic, ::Val{true})
    normalform!(α)
    normalform!(β)
    return coeffs(α) == coeffs(β)
end

function Base.:(==)(α::Cyclotomic, β::Cyclotomic, ::Val{false})
    l = lcm(conductor(α), conductor(β))
    return ==(embed(α, l), embed(β, l), Val{true}())
end

Base.iszero(α::Cyclotomic) =
    all(iszero, values(α)) || (normalform!(α); all(iszero, values(α)))
Base.isone(α::Cyclotomic) = throw("Not implemented")

function isnormalized(α::Cyclotomic, basis)
    # return all(in(basis), exponents(α))
    for e in exponents(α)
        !(e in basis) && return false
    end
    return true
end

function _all_equal(α::Cyclotomic, exps, value)
    # return all(==(value), (α[e] for e in exps))
    for e in exps
        α[e] == value || return false
    end
    return true
end
