####
#   Arithmetic

zero!(α::Cyclotomic{T}) where T = (coeffs(α) .= zero(T); α)
one!(α::Cyclotomic{T}) where T = (zero!(α); α[0] = one(α[0]); α)
Base.zero(α::Cyclotomic, m::Integer=conductor(α)) = zero!(similar(α, m))
Base.one(α::Cyclotomic) = one!(similar(α))

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
    copyto!(coeffs(out), coeffs(mul!(dense(out), α, β)))
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

    for (αe, αc) in α
        for (βe, βc) in β
            out[αe+βe] += αc * βc
        end
    end

    return out
end

for (op, fn) in ((:+, :add!), (:-, :sub!), (:*, :mul!))
    @eval begin
        function Base.$op(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
            if conductor(α)==conductor(β)
                return $fn(similar(α, promote_type(T, S)), α, β)
            else
                l = lcm(conductor(α), conductor(β))
                return $op(embed(α, l), embed(β, l))
            end
        end
    end
end

function Base.conj!(out::Cyclotomic, α::Cyclotomic, n::Integer=-1)
    zero!(out)
    for (exp, c) in α
        out[n*exp] = c
    end
    return out
end

function Base.conj(α::Cyclotomic, n::Integer=-1)
    return conj!(similar(α), α, n)
end

galois_conj(α::Cyclotomic, n::Integer=-1) =
    (@assert gcd(n, conductor(α)) == 1; conj(α, n))

function inv!(out::Cyclotomic{T}, α::Cyclotomic) where T
    copyto!(coeffs(out), coeffs(inv!(dense(out), α)))
    return out
end

function inv!(out::Cyclotomic{T, <:DenseVector}, α::Cyclotomic, tmp=similar(out)) where T
    if out === α
        out = similar(out)
    end

    one!(out)
    tmp2 = deepcopy(out)

    basis, fb = zumbroich_viacomplement(conductor(α))
    lb = length(basis)
    conjugates_counter = 0

    for i in 2:conductor(α)-1
        conjugates_counter == lb-1 && break # break finish
        all(x->gcd(i, first(x)) == 1, fb) || continue
        conjugates_counter += 1
        mul!(tmp2, out, conj!(tmp, α, i))
        copyto!(coeffs(out), coeffs(tmp2))
    end

    # out is now the product of non-trivial Galois conjugates of α:
    # out = Π_{σ(Gal(ℚ(ζ_n)/ℚ)), σ≠id} σ(α)
    # since Π_{σ(Gal(ℚ(ζ_n)/ℚ))} σ(α) = norm_ℚ(α) ∈ ℚ we have
    # 1 = α·out/(α·out) = α · out/norm_ℚ(α), hence
    # α¯¹ = out/norm_ℚ(α)

    norm_ℚ = reduced_embedding(mul!(tmp2, out,α))
    @assert conductor(norm_ℚ) == 1 # norm_ℚ is real
    norm_α = norm_ℚ[0]

    out = mul!(out, out, inv(norm_α))

    return out
end

Base.inv(α::Cyclotomic) = inv!(similar(α), α)
