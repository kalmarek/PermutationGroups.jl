####
#   Arithmetic

zero!(α::Cyclotomic{T}) where T = (coeffs(α) .= zero(T); α)
Base.zero(α::Cyclotomic, m::Integer=conductor(α)) =
    (res = similar(α, m); zero!(res))
Base.one(α::Cyclotomic{T}) where T = (res = zero(α); res[0] = one(T); res)

############################
# Module structure:

Base.:-(α::Cyclotomic) = Cyclotomic(-coeffs(α))

for op in (:+, :-)
    @eval begin
        function Base.$op(α::Cyclotomic{T}, r::R) where {T, R<:Real}
            res = similar(α, promote_type(T, R))
            copyto!(coeffs(res), coeffs(α))
            res[0] = $op(res[0], r)
            return res
        end
    end
end

Base.:+(r::Real, α::Cyclotomic) = α + r
Base.:-(r::Real, α::Cyclotomic) = (res = -α; res[0] += r; res)

mul!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= coeffs(α) .* c; out)
div!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= div.(coeffs(α), c); out)

Base.:*(c::T, α::Cyclotomic{S}) where {S,T<:Real} =
    mul!(similar(α, promote_type(S,T)), α, c)
Base.:*(α::Cyclotomic, c::T) where T<:Real = c*α
Base.:(//)(α::Cyclotomic, c::Real) = Cyclotomic(coeffs(α).//c)
Base.:(/)(α::Cyclotomic, c::Real) = Cyclotomic(coeffs(α)./c)

function Base.div(α::Cyclotomic, c::Number)
    res = similar(α, first(Base.return_types(div, (valtype(α), typeof(c)))))
    return div!(res, α, c)
end

###########################
# Ring structure:

add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .+ coeffs(β); out)
sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .- coeffs(β); out)

function mul!(out::Cyclotomic{T}, α::Cyclotomic, β::Cyclotomic) where T
    tmp = Cyclotomic{T, Vector{T}}(conductor(out), zeros(T, conductor(out)))
    mul!(tmp, α, β)
    copyto!(coeffs(out), coeffs(tmp))
    return out
end

function mul!(
    out::Cyclotomic{T,<:DenseVector},
    α::Cyclotomic,
    β::Cyclotomic,
) where {T}
    if out === α || out === β
        out = similar(out)
    end
    zero!(out)

    for (αe, αc) in zip(exponents(α), values(α))
        for (βe, βc) in zip(exponents(β), values(β))
            out[αe+βe] += αc * βc
        end
    end

    return out
end

for (op, fn) in ((:+, :add!), (:-, :sub!), (:*, :mul!))
    @eval begin
        Base.$op(α::Cyclotomic, β::Cyclotomic) =
            $op(α, β, Val{conductor(α)==conductor(β)}())

        function Base.$op(α::Cyclotomic{T}, β::Cyclotomic{S}, ::Val{true}) where {T,S}
            return $fn(similar(α, promote_type(T, S)), α, β)
        end

        function Base.$op(α::Cyclotomic, β::Cyclotomic, ::Val{false})
            l = lcm(conductor(α), conductor(β))
            return $op(embed(α, l), embed(β, l), Val{true}())
        end
    end
end

function Base.conj!(out::Cyclotomic, α::Cyclotomic, n::Integer=-1)
    zero!(out)
    for (exp, c) in zip(exponents(α), values(α))
        out[n*exp] = c
    end
    return out
end

function Base.conj(α::Cyclotomic, n::Integer=-1)
    return conj!(similar(α), α, n)
end

galois_conj(α::Cyclotomic, n::Integer=-1) =
    (@assert gcd(n, conductor(α)) == 1; conj(α, n))

