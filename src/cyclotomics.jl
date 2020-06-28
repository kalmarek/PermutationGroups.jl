module Cyclotomics

using Primes
using SparseArrays

export zumbroich, Cyclotomic, exponents, conductor, E


###############################################################################
#
#   Zumproich basis implementation

function J(k::Integer, p::Integer)
    if p == 2
        return ifelse(iszero(k), 0:0, 0:1)
    else
        pm1 = convert(Int, p) - 1
        return ifelse(iszero(k), 1:pm1, -pm1>>1:pm1>>1)
    end
end

function zumbroich_plain(n::Integer, m::Integer = 1)

    # following the description of
    # https://www.gap-system.org/Manuals/doc/ref/chap60_mj.html#X7F52BEA0862E06F2

    @assert iszero(last(divrem(n, m)))
    isone(n) && return [0]

    factors_n = factor(n)
    factors_m = factor(m)

    exps = Vector{Int}[]
    sizehint!(exps, 2 * length(factors_n))

    for (p, ν) in factors_n
        for k = factors_m[p]:ν-1
            push!(exps, [div(n, p^(k + 1)) * j for j in J(k, p)])
        end
    end

    w = (sum.(Iterators.product(exps...))) .% n
    return sort!(vec(w))
end

###############################################################################
#
# Forbidden residues, i.e. the complement of the Zumbroich basis

function _forbidden_residues(n, p, ν, q = p^ν)
    k = div(n, q)
    if isodd(p)
        a = (p^(ν - 1) - 1) >> 1
        return BitSet((k*(-a+q):k:k*(a+q)) .% q)
    else # iseven(p)
        return BitSet((k*(q>>1):k:k*(q-1)) .% q)
    end
end

struct ForbiddenResidues{I}
    primes_powers::Vector{Tuple{I,I,BitSet}}
end

function ForbiddenResidues(
    n::Integer,
    pf::Primes.Factorization{I} = factor(n),
) where {I}
    l = length(pf)
    primes_powers = Vector{Tuple{I,I,BitSet}}(undef, l)
    for (i, (p, ν)) in enumerate(pf)
        pν = p^ν
        primes_powers[i] = (p, pν, _forbidden_residues(n, p, ν, pν))
    end
    return ForbiddenResidues(primes_powers)
end

Base.length(fr::ForbiddenResidues) = length(fr.primes_powers)

Base.iterate(fr::ForbiddenResidues, s = 1) =
    (s > length(fr) ? nothing : (fr.primes_powers[s], s + 1))

@inline function Base.in(n::Integer, fr::ForbiddenResidues)
    @assert n ≥ 0
    for (_, q, residues) in fr
        n % q ∈ residues && return true
    end
    return false
end

function _zumbroich_complement_old(n::Integer, factor_n = factor(n))
    # following the wonderful documenting comments at the top of
    # https://github.com/gap-system/gap/blob/master/src/cyclotom.c

    factor_n = factor(n)
    forbidden =
        [(p, p^ν) => _forbidden_residues(n, p, ν, p^ν) for (p, ν) in factor_n]
    isone(n) && return [0], forbidden

    exps = Vector{typeof(n)}(undef, totient(factor_n))
    count = 0
    for i = 0:n-1
        for ((_, q), forbidden_residues) in forbidden
            if i % q ∈ forbidden_residues
                @goto next_i
            end
        end
        count += 1
        exps[count] = i
        @label next_i
    end
    @assert count == length(exps)
    return exps, forbidden
end


function _zumbroich_complement(n::Integer, factor_n = factor(n))

    # following the wonderful documenting comments at the top of
    # https://github.com/gap-system/gap/blob/master/src/cyclotom.c

    forbidden = ForbiddenResidues(n, factor_n)

    isone(n) && return [0], forbidden

    exps = Vector{typeof(n)}(undef, totient(factor_n))
    count = 0
    for i = 0:n-1
        i ∈ forbidden && continue
        count += 1
        exps[count] = i
        # count == length(exps) && break
    end
    @assert count == length(exps)
    return exps, forbidden
end

zumbroich_complement(n::Integer) = first(_zumbroich_complement(n))

function zumbroich_direct(n::Integer)
    isone(n) && return [0]
    basis = [0]

    factor_n = factor(n)

    if iseven(n)
        for k = 1:factor_n[2]-1
            basis = union!(basis, basis .+ (n >> (k + 1)))
        end
    end

    for (p, ν) in factor_n
        p == 2 && continue
        for k = 0:ν-1
            ran = div(n, p^(k + 1)) .* J(k, p)
            new_basis = typeof(basis)()
            sizehint!(new_basis, length(basis) * length(ran))
            for b in basis
                append!(new_basis, ran .+ b)
            end
            basis = new_basis
        end
    end

    for (i, x) in enumerate(basis)
        basis[i] = mod(basis[i], n)
    end
    unique!(basis)
    sort!(basis)
    return basis
end

zumbroich = zumbroich_complement

###############################################################################
#
#   Cyclotomics
"""
    Cyclotomic(n, coeffs::AbstractVector)
Element of `n`-th cyclotomic field with coefficients stored as `coeffs`.

To access the internals of a cyclotomic use API functions:
 * `conductor` - the conductor of a cyclotomic, i.e. the `n` used currently for storage. This might not be the minimal embeding field of a cyclotomic.
 * `getindex`/`setindex!` - use `α[i]` to access the coefficient at `i`-th power of a cyclotomic (in a circular fashion)
 * `values`/`exponents` - paired iterators over _non zero_ coefficients/exponents corresponding to _non-zero_ coefficients
 * `normalform!` - bring a cyclotomic to its unique representation as given by Zumbroich basis (also available in non-modifying form).

!!! warning "Beware!"

    `hash` function will not reduce a cyclotomic to its minimal embedding field, as this may be a very expensive operation, and will compute the `hash` of a cyclotomic _at current embeding_. Therefore _equal cyclotomics_ in different embeddings may have _different hashes_! To avoid this pitfall use `normalform!(α, minimalfield(α))`.
"""
struct Cyclotomic{T, A<:AbstractVector{T}} <: Number
    n::Int
    coeffs::A
end

Cyclotomic(v::V) where V<:AbstractVector = Cyclotomic{eltype(v), V}(length(v), v)
Cyclotomic{T}(α::Cyclotomic) where T = Cyclotomic(conductor(α), convert.(T, α.coeffs))

"""
    E(n[, i=1])
Return the `i`-th power of `n`-th root of unity with sparse vector as storage.
"""
function E(n, i=1)
    k = totient(n)
    i = (0 <= i < n ? i : mod(i, n))
    coeffs = sparsevec([i+1], [1], n);
    sizehint!(coeffs.nzind, k)
    sizehint!(coeffs.nzval, k)
    return Cyclotomic(n, coeffs)
end

####
#   Low-level interface

conductor(α::Cyclotomic) = α.n
coeffs(α::Cyclotomic) = α.coeffs

_to_index(α::Cyclotomic, idx::Integer) = mod(idx, conductor(α)) + 1

Base.@propagate_inbounds function Base.getindex(α::Cyclotomic, exp::Integer)
    return coeffs(α)[_to_index(α, exp)]
end

Base.getindex(α::Cyclotomic, itr) = [α[i] for i in itr]

Base.@propagate_inbounds function Base.setindex!(α::Cyclotomic, val, exp::Integer)
    coeffs(α)[_to_index(α, exp)] = val
    return val
end

Base.firstindex(::Cyclotomic) = 0
Base.lastindex(α::Cyclotomic) = conductor(α) - 1
Base.eachindex(α::Cyclotomic) = firstindex(α):lastindex(α)

# general definitions for iteration
Base.values(α::Cyclotomic) = (c for c in coeffs(α) if !iszero(c))
exponents(α::Cyclotomic) = (e for e in eachindex(α) if !iszero(α[e]))

# sparse storage definitions
Base.values(α::Cyclotomic{T, <:SparseVector}) where T =
    (c for c in coeffs(α).nzval if !iszero(c))
exponents(α::Cyclotomic{T, <:SparseVector}) where T =
    (e-1 for e in coeffs(α).nzind if !iszero(α[e-1]))

Base.valtype(::Type{Cyclotomic{T}}) where T = T
Base.valtype(::Cyclotomic{T}) where T = T

Base.similar(α::Cyclotomic, T=valtype(α)) = Cyclotomic(similar(coeffs(α), T))

Base.deepcopy_internal(α::Cyclotomic, ::IdDict) = Cyclotomic(deepcopy(coeffs(α)))

function Base.hash(α::Cyclotomic, h::UInt)
    normalform!(α)
    return hash(coeffs(α), hash(conductor(α), hash(Cyclotomic, h)))
end

function Base.:(==)(α::Cyclotomic, β::Cyclotomic)
    conductor(α) == conductor(β) || throw("Embeddings are not implemented yet")
    normalform!(α)
    normalform!(β)
    return coeffs(α) == coeffs(β)
end

Base.iszero(α::Cyclotomic) = all(iszero, values(α)) || (normalform!(α); all(iszero, values(α)))
Base.isone(α::Cyclotomic) = throw("Not implemented")

####
#   normalform (reduction to Zumbroich basis)

normalform(α::Cyclotomic) = normalform!(deepcopy(α))

function normalform!(α::Cyclotomic{T}, tmp=Cyclotomic{T, Vector{T}}(conductor(α), coeffs(α))) where T
    n = conductor(α)

    basis, forbidden = let (b, fb) =_zumbroich_complement(n)
        BitSet(b), fb
    end

    all(in(basis), exponents(α)) && return α

    for fb in forbidden
        for exp in exponents(tmp)
            exp in basis && continue
            _replace_exponent!(tmp, exp, fb)
        end
    end

    α.coeffs .= coeffs(tmp)
    return α
end

function _replace_exponent!(α::Cyclotomic, exp::Integer, fb)
    val = α[exp]
    p, q, forbidden_exps = fb

    # @debug p forbidden_exps
    if exp % q ∈ forbidden_exps
        # @debug "removing $(α[exp])*ζ^$exp from" α
        m = conductor(α) ÷ p
        for i in 1:p-1
            α[exp + i*m] -= val
            # exp + i*m is not guaranteed to be allowed, but it's larger than exp and hence will be dealt later
        end
        α[exp] = 0
        # @debug "after removal:" α
    end
    return α
end

####
#   Arithmetic

zero!(α::Cyclotomic) = (coeffs(α) .= 0; α)
Base.zero(α::Cyclotomic) = (res = similar(α); zero!(res))
Base.one(α::Cyclotomic) = (o = zero(α); o[0] = 1; o)

add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .+ coeffs(β); out)
sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (coeffs(out) .= coeffs(α) .- coeffs(β); out)
mul!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= coeffs(α) .* c; out)
div!(out::Cyclotomic, α::Cyclotomic, c::Real) =
    (coeffs(out) .= div.(coeffs(α), c); out)

add!(α::Cyclotomic, β::Cyclotomic) = add!(α, α, β)
sub!(α::Cyclotomic, β::Cyclotomic) = sub!(α, α, β)
mul!(α::Cyclotomic{T}, c::T) where T = mul!(α, α, c)
div!(α::Cyclotomic{T}, c::T) where T = div!(α, α, c)

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

####
#   IO

function subscriptify(n::Int)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join((Char(subscript_0 + i) for i in reverse(digits(n))))
end

function superscriptify(n::Int)
    superscripts = Dict(
    0 => "⁰",
    1 => "¹",
    2 => "²",
    3 => "³",
    4 => "⁴",
    5 => "⁵",
    6 => "⁶",
    7 => "⁷",
    8 => "⁸",
    9 => "⁹",
)
    return join(superscripts[d] for d in reverse(digits(n)))
end

function Base.show(io::IO, α::Cyclotomic{T}) where T
    ζ = "ζ"*subscriptify(conductor(α))
    if iszero(α)
        print(io, zero(T))
    else
        for (i, exp) in enumerate(exponents(α))
            val = α[exp]
            if iszero(exp)
                print(io, val)
                continue
            end
            sign_str = (val >= zero(T) ? (i == 1 ? "" : " +") : " ")
            exp_str = isone(exp) ? "" : "$(superscriptify(exp))"
            print(io, sign_str, val, "*", ζ, exp_str)
        end
    end
end

end # of module Cyclotomics
