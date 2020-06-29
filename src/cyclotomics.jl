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

function _to_index(α::Cyclotomic, idx::Integer)
    # return mod(idx, conductor(α)) + 1
    0 <= idx < conductor(α) && return idx + 1
    conductor(α) <= idx && return (idx % conductor(α))+1
    return (idx % conductor(α)) + conductor(α) + 1
end

Base.@propagate_inbounds function Base.getindex(α::Cyclotomic, exp::Integer)
    return α.coeffs[_to_index(α, exp)]
end

Base.getindex(α::Cyclotomic, itr) = [α[i] for i in itr]

Base.@propagate_inbounds function Base.setindex!(α::Cyclotomic, val, exp::Integer)
    α.coeffs[_to_index(α, exp)] = val
    return val
end

Base.@propagate_inbounds function Base.setindex!(α::Cyclotomic, val, itr)
    for exp in itr
        α[exp] = val
    end
    return itr
end

Base.firstindex(::Cyclotomic) = 0
Base.lastindex(α::Cyclotomic) = conductor(α) - 1
Base.eachindex(α::Cyclotomic) = firstindex(α):lastindex(α)

# general definitions for iteration
Base.values(α::Cyclotomic) = (α[i] for i in eachindex(α) if !iszero(α[i]))
exponents(α::Cyclotomic)   = (i    for i in eachindex(α) if !iszero(α[i]))

# sparse storage definitions
Base.values(α::Cyclotomic{T, <:SparseVector}) where T =
    (c for c in coeffs(α).nzval if !iszero(c))
exponents(α::Cyclotomic{T, <:SparseVector}) where T =
    (e-1 for (i,e) in enumerate(coeffs(α).nzind) if !iszero(coeffs(α).nzval[i]))

Base.valtype(::Type{Cyclotomic{T}}) where T = T
Base.valtype(::Cyclotomic{T}) where T = T

Base.similar(α::Cyclotomic, T::Type=valtype(α)) = similar(α, T, conductor(α))
Base.similar(α::Cyclotomic, m::Integer) = similar(α, valtype(α), m)
Base.similar(α::Cyclotomic, T::Type, n::Integer) = Cyclotomic(similar(coeffs(α), T, n))

function Base.hash(α::Cyclotomic, h::UInt)
    normalform!(α)
    return hash(coeffs(α), hash(conductor(α), hash(Cyclotomic, h)))
end

Base.:(==)(α::Cyclotomic, β::Cyclotomic) =
    ==(α, β, Val{conductor(α)==conductor(β)}())

function Base.:(==)(α::Cyclotomic, β::Cyclotomic, ::Val{true})
    normalform!(α)
    normalform!(β)
    return coeffs(α) == coeffs(β)
end

function Base.:(==)(α::Cyclotomic, β::Cyclotomic, ::Val{false})
    l = lcm(conductor(α), conductor(β))
    return ==(embed(α, l), embed(β, l), Val{true}())
end

Base.iszero(α::Cyclotomic) = all(iszero, values(α)) || (normalform!(α); all(iszero, values(α)))
Base.isone(α::Cyclotomic) = throw("Not implemented")

####
#   normalform (reduction to Zumbroich basis)

normalform(α::Cyclotomic) = normalform!(deepcopy(α))

function isnormalized(α::Cyclotomic, basis)
    # return all(in(basis), exponents(α))
    for e in exponents(α)
        !(e in basis) && return false
    end
    return true
end

function normalform!(α::Cyclotomic{T};
    tmp=Cyclotomic{T, Vector{T}}(conductor(α), coeffs(α)),
    basis_forbidden = _zumbroich_complement(conductor(α))
    ) where T
    # @debug "normalizing:" conductor(α) coeffs(α)

    basis, forbidden = basis_forbidden

    isnormalized(α, BitSet(basis)) && return α

    tmp.coeffs .= coeffs(α)

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
        # @debug " " collect(α[exp .+ 1:m:m*(p-1)])

        # TODO: this doesn't work
        # α[exp .+ 1:m:m*(p-1)] .-= val

        for i in 1:p-1
            α[exp + i*m] -= val
            # exp + i*m is not guaranteed to be allowed, but it's larger than exp and hence will be dealt later
        end
        α[exp] = 0
        # @debug "after removal:" α
    end
    return α
end

function embed(α::Cyclotomic, m::Integer)
    conductor(α) == m && return deepcopy(α)
    if conductor(α) > m
        @assert conductor(α) % m == 0 "Embeding of ℚ(ζ$(subscriptify(m))) ↪ ℚ(ζ$(subscriptify(conductor(α)))) is not possible."
        return reduced_embedding(α, m)
    else
        k, _r = divrem(m, conductor(α))
        @assert _r == 0 "Embeding of ℚ(ζ$(subscriptify(conductor(α)))) ↪ ℚ(ζ$(subscriptify(m))) is not possible."

        res = zero!(similar(α, valtype(α), m))

        @inbounds for e in exponents(α)
            res[e*k] = α[e]
        end
        return res
    end
end

function _all_equal(α::Cyclotomic, exps, value)
    # return all(==(value), (α[e] for e in exps))
    for e in exps
        α[e] == value || return false
    end
    return true
end

function reduced_embedding(α::Cyclotomic{T,V}, m::Integer=1) where {T, V}
    k = gcd(exponents(α)...)
    @debug "gcd(exponents(α)) = $k" collect(exponents(α))

    tmp = let (d,r) = divrem(conductor(α), k)
        if (k > 1 && iszero(r))
            @debug "Reducing the embedding ring to $T(ζ$(subscriptify(d)))"
            tmp = Cyclotomic{T, Vector{T}}(d, zeros(T, d))
        else
            n = conductor(α)
            tmp = Cyclotomic{T, Vector{T}}(n, zeros(T, n))
        end
        tmp # tmp has dense storage
    end

    if k > 1
        # the trivial reduction E(n)^e → E(n÷k)^(e÷k)
        @debug "Performing trivial reduction from ℚ(ζ$(subscriptify(conductor(α)))) → ℚ(ζ$(subscriptify(conductor(tmp))))"
        @inbounds for (e,v) in zip(exponents(α), values(α))
            tmp[div(e,k)] = α[e]
        end
    else
        @debug "No trivial reduction is possible"
        tmp.coeffs .= coeffs(α)
    end

    normalform!(tmp)

    basis, forbidden = Cyclotomics._zumbroich_complement(conductor(tmp))

    phi_nc = length(basis)
    nz = count(!iszero, coeffs(tmp))
    val = first(values(tmp))

    if all(p == q for (p,q,_) in forbidden) && # conductor(tmp) is square-free
        phi_nc == nz && # tmp is supported on φ(n) elements
        _all_equal(tmp, exponents(tmp), val) # all nz coeffs of tmp are equal

        @debug "Cyclotomic is real:" α = -val
        res = Cyclotomic{T, V}(m, similar(coeffs(α), m))
        zero!(res)
        res[0] = -val
        return res
    end

    n = conductor(tmp)
    for (p,q,fb_res) in forbidden

        p != q && continue
        @debug "$(p)² does NOT divide conductor $(conductor(tmp))" p q

        m % p == 0 && continue
        @debug "prime $p does not divide $m so embedding is possible"

        # there are nc/p residue classes, each contains p-1 elts, and tmp is constant on each
        nnz, _r = divrem(nz, p-1)
        @debug "the number of non-zero coeffs $nz ≡ $_r mod $(p-1)"
        _r != 0 && continue


        n╱p = n ÷ p

        # check that coeffs(tmp) are constant on each congruence classes
        equal_on_classes = true
        for i in 1:p:n
            if !_all_equal(tmp, i.+(2n╱p:n╱p:n+i), tmp[i + n╱p])
                equal_on_classes = false
                @debug "skipping reduction over prime $p as cyclotomic is not constant over class $(i+(n╱p)) (mod $(n╱p)) " tmp[i.+(2n╱p:n╱p:n+i)]
                break
            end
        end
        equal_on_classes || continue # to next prime

        @debug "cyclotomic is constant over residue classes (mod $(n╱p)), reducing"


        for i in 1:p:n # inverse to normalform over p
            tmp[i] = -tmp[i+n╱p]
            tmp[(i+n╱p):n╱p:(n╱p+i)] = zero(T)
        end

        for i in 1:n╱p
            tmp[i] = tmp[i*p]
            tmp[i*p] = zero(T)
        end

        nz = nnz # update the number on non-zero coefficients
        tmp = Cyclotomic{T, Vector{T}}(np, resize!(coeffs(tmp), np))
    end


    return Cyclotomic{T, V}(conductor(tmp), coeffs(tmp))
end

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
            val_str = "$val"
            exp_str = isone(exp) ? "" : "$(superscriptify(exp))"
            print(io, sign_str, val_str, "*", ζ, exp_str)
        end
    end
end

end # of module Cyclotomics
