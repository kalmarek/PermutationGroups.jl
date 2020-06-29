####
#   Arithmetic

zero!(α::Cyclotomic{T}) where T = (coeffs(α) .= zero(T); α)
Base.zero(α::Cyclotomic, m::Integer=conductor(α)) =
    (res = similar(α, m); zero!(res))
Base.one(α::Cyclotomic{T}) where T = (res = zero(α); res[0] = one(T); res)

add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .+ coeffs(β); out)
sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .- coeffs(β); out)
mul!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= coeffs(α) .* c; out)
div!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= div.(coeffs(α), c); out)

function Base.:+(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
    conductor(α) == conductor(β) || throw("Not Implemented")
    return add!(similar(α, promote_type(T, S)), α, β)
end

function Base.:-(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
    conductor(α) == conductor(β) || throw("Not Implemented")
    return sub!(similar(α, promote_type(T, S)), α, β)
end

for op in (:+, :-)
    @eval begin
        function Base.$op(α::Cyclotomic{T}, r::R) where {T, R<:Real}
            res = similar(α, promote_type(T, R))
            coeffs(res) .= coeffs(α)
            res[0] = $op(res[0], r)
            return res
        end
    end
end

Base.:+(r::Real, α::Cyclotomic) = α + r
Base.:-(r::Real, α::Cyclotomic) = (res = -α; res[0] += r; res)
Base.:-(α::Cyclotomic) = Cyclotomic(-coeffs(α))
Base.:*(c::T, α::Cyclotomic{S}) where {S,T<:Number} =
    mul!(similar(α, promote_type(S,T)), α, c)
Base.:*(α::Cyclotomic, c::T) where T<:Number = c*α
Base.:(//)(α::Cyclotomic, c::Number) = Cyclotomic(coeffs(α).//c)
Base.:(/)(α::Cyclotomic, c::Number) = Cyclotomic(coeffs(α)./c)

function Base.div(α::Cyclotomic, c::Number)
    res = similar(α, first(Base.return_types(div, (valtype(α), typeof(c)))))
    return div!(res, α, c)
end

function mul!(out::Cyclotomic{T}, α::Cyclotomic, β::Cyclotomic) where T
    if out === α || out === β
        out = deepcopy(out)
    end
    tmp = Cyclotomic{T, Vector{T}}(conductor(out), zeros(T, conductor(out)))

    for (αe, αc) in zip(exponents(α), values(α))
        for (βe, βc) in zip(exponents(β), values(β))
            tmp[αe + βe] += αc*βc
        end
    end
    zero!(out)
    out.coeffs .= coeffs(tmp)
    return out
end

Base.:*(α::Cyclotomic, β::Cyclotomic) = mul!(similar(α), α, β)
