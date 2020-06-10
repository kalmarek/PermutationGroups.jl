module Cyclotomics

using Primes

export zumbroich, Cyclotomic, exponents, conductor


###############################################################################
#
#   Zumproich basis implementation

function J(k::Integer, p::Integer)
    if p == 2
        return ifelse(iszero(k), 0:0, 0:1)
    else
        pm1 = convert(Int, p) - 1
        return ifelse(iszero(k), 1:pm1, -pm1>>1:pm1>>1 )
    end
end

function zumbroich_plain(n::Integer, m::Integer=1)

    # following the description of
    # https://www.gap-system.org/Manuals/doc/ref/chap60_mj.html#X7F52BEA0862E06F2

    @assert iszero(last(divrem(n, m)))
    isone(n) && return [0]

    factors_n = factor(n)
    factors_m = factor(m)

    exps = Vector{Int}[]
    sizehint!(exps, 2length(factors_n))

    for (p, ν) in factors_n
        for k in factors_m[p]:ν-1
            push!(exps, [div(n, p^(k+1))*j for j in J(k,p)])
        end
    end

    w = (sum.(Iterators.product(exps...))) .% n
    return sort!(vec(w))
end

###############################################################################
#
# Forbidden residues, i.e. the complement of the Zumbroich basis

function _forbidden_residues(n, p, ν, q=p^ν)
    k = div(n, q)
    if isodd(p)
        a = (p^(ν-1) - 1) >> 1
        return BitSet((k*(-a+q):k:k*(a+q)) .% q)
    else # iseven(p)
        return BitSet((k*(q >> 1):k:k*(q-1)) .% q)
    end
end

struct ForbiddenResidues{I}
    primes_powers::Vector{Tuple{I,I,BitSet}}
end

function ForbiddenResidues(n::Integer, pf::Primes.Factorization{I}=factor(n)) where I
    l = length(pf)
    primes_powers = Vector{Tuple{I,I, BitSet}}(undef, l)
    for (i, (p,ν)) in enumerate(pf)
        pν = p^ν
        primes_powers[i] = (p, pν, _forbidden_residues(n, p, ν, pν))
    end
    return ForbiddenResidues(primes_powers)
end

Base.length(fr::ForbiddenResidues) = length(fr.primes_powers)

Base.iterate(fr::ForbiddenResidues, s=1) =
    (s > length(fr) ? nothing : (fr.primes_powers[s], s+1))

@inline function Base.in(n::Integer, fr::ForbiddenResidues)
    @assert n ≥ 0
    for (_, q, residues) in fr
        n % q ∈ residues && return true
    end
    return false
end

function _zumbroich_complement_old(n::Integer, factor_n=factor(n))
    # following the wonderful documenting comments at the top of
    # https://github.com/gap-system/gap/blob/master/src/cyclotom.c

    factor_n = factor(n)
    forbidden = [(p,p^ν)=>_forbidden_residues(n, p, ν, p^ν) for (p, ν) in factor_n]
    isone(n) && return [0], forbidden

    exps = Vector{typeof(n)}(undef, totient(factor_n))
    # @show forbidden
    count = 0
    for i in 0:n-1
        for ((_,q), forbidden_residues) in forbidden
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


function _zumbroich_complement(n::Integer, factor_n=factor(n))

    # following the wonderful documenting comments at the top of
    # https://github.com/gap-system/gap/blob/master/src/cyclotom.c

    forbidden = ForbiddenResidues(n, factor_n)

    isone(n) && return [0], forbidden

    exps = Vector{typeof(n)}(undef, totient(factor_n))
    count = 0
    for i in 0:n-1
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
        for k in 1:factor_n[2]-1
            basis = union!(basis, basis .+ (n>>(k+1)))
        end
    end

    for (p, ν) in factor_n
        p == 2 && continue
        for k in 0:ν-1
            ran = div(n, p^(k+1)).*J(k, p)
            new_basis = typeof(basis)()
            sizehint!(new_basis, length(basis)*length(ran))
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

struct Cyclotomic{T, A<:AbstractVector{T}} <: Number
    coeffs::A
end

Cyclotomic{T}(α::Cyclotomic) where T = Cyclotomic(convert.(T, α.coeffs))
Cyclotomic(v::V) where V<:AbstractVector = Cyclotomic{eltype(v), V}(v)

function E(i, n)
    k = totient(n)
    i = (0 <= i < n ? i : mod(i, n))
    coeffs = sparsevec([i+1], [1], n);
    sizehint!(coeffs.nzind, k)
    sizehint!(coeffs.nzval, k)
    return Cyclotomic(coeffs)
end

####
#   Low-level interface

conductor(α::Cyclotomic) = length(α.coeffs)

_to_index(α::Cyclotomic, idx::Integer) = mod(idx, conductor(α)) + 1

Base.@propagate_inbounds function Base.getindex(α::Cyclotomic, exp::Integer)
    return α.coeffs[_to_index(α, exp)]
end

Base.getindex(α::Cyclotomic, itr) = [α[i] for i in itr]

Base.@propagate_inbounds function Base.setindex!(α::Cyclotomic, val, exp::Integer)
    α.coeffs[_to_index(α, exp)] = val
    return val
end

# general definitions
Base.values(α::Cyclotomic) = (c for c in α.coeffs if !iszero(c))
exponents(α::Cyclotomic) = (e for e in 0:conductor(α)-1 if !iszero(α[e]))

# sparse storage definitions
Base.values(α::Cyclotomic{T, <:SparseVector}) where T =
    (c for c in α.coeffs.nzval if !iszero(c))
exponents(α::Cyclotomic{T, <:SparseVector}) where T =
    (e-1 for (i,e) in enumerate(α.coeffs.nzind) if !iszero(α[e-1]))

Base.similar(α::Cyclotomic, T::DataType=eltype(α)) = Cyclotomic(similar(α.coeffs, T))

Base.eltype(::Type{Cyclotomic{T}}) where T = T
Base.eltype(α::Cyclotomic{T}) where T = T

Base.deepcopy_internal(α::Cyclotomic, ::IdDict) = Cyclotomic(deepcopy(α.coeffs))

# TODO: hash, predicates: ==, iszero, isone
Base.iszero(α::Cyclotomic) = all(iszero, values(α))
# Base.isone(α::Cyclotomic) = isone(c[0]) && all(iszero(α.coeffs[i]) for i in 1:length(conductor))

####
#   normalform (reduction to Zumbroich basis)

normalform(α::Cyclotomic, n=conductor(α)) = normalform!(deepcopy(α), n)

function normalform!(α::Cyclotomic, n=conductor(c))
    n == conductor(α) || throw("Embeddings are not implemented yet")

    basis, forbidden = let
        b, f = Cyclotomics._zumbroich_complement(n)
        BitSet(b), f
    end
    @show forbidden

    all(in(basis), exponents(α)) && return α

    for ((p, q), forbidden_residues) in forbidden
        # @debug p forbidden_residues
        for i in exponents(α)
            if i % q ∈ forbidden_residues
                # @debug "removing $(α[i])*ζ^$i from" α
                ndivp = div(n, p)
                for e in 1:p-1
                    α[i+e*ndivp] -= α[i]
                end
                α[i] = 0
                # @debug "after removal:" α
            end
        end
    end
    return α
end

####
#   Arithmetic

Base.zero!(α::Cyclotomic) = (α.coeffs .= 0; α)
Base.zero(α::Cyclotomic) = (res = similar(α); zero!(res))
Base.one(α::Cyclotomic) = (o = zero(α); o[0] = 1; o)

AbstractAlgebra.add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (out.coeffs .= α.coeffs .+ β.coeffs; out)
sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic) =
    (out.coeffs .= α.coeffs .- β.coeffs; out)
AbstractAlgebra.mul!(out::Cyclotomic, α::Cyclotomic, c::Number) =
    (out.coeffs .= α.coeffs .* c; out)
div!(out::Cyclotomic, α::Cyclotomic, c::Number) =
    (out.coeffs .= div.(α.coeffs, c); out)

AbstractAlgebra.add!(α::Cyclotomic, β::Cyclotomic) = add!(α, α, β)
sub!(α::Cyclotomic, β::Cyclotomic) = sub!(α, α, β)
AbstractAlgebra.mul!(α::Cyclotomic{T}, c::T) where T = mul!(α, α, c)
div!(α::Cyclotomic{T}, c::T) where T = div!(α, α, c)

function Base.:+(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
    conductor(α) == conductor(β) || throw("Not Implemented")
    return add!(similar(α, promote_type(T, S)), α, β)
end

function Base.:-(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
    conductor(α) == conductor(β) || throw("Not Implemented")
    return sub!(similar(α, promote_type(T, S)), α, β)
end

function Base.:+(α::Cyclotomic{T}, r::R) where {T, R<:Real}
    res = similar(α, promote_type(T, R))
    res.coeffs .= α.coeffs
    res.coeffs[0] += r
    return res
end

function Base.:-(α::Cyclotomic{T}, r::R) where {T, R<:Real}
    res = similar(α, promote_type(T, R))
    res.coeffs .= α.coeffs
    res.coeffs[0] -= r
    return res
end

Base.:-(α::Cyclotomic) = Cyclotomic(-α.coeffs)
Base.:*(c::T, α::Cyclotomic{S}) where {S,T<:Number} =
    mul!(similar(α, promote_type(S,T)), α, c)
Base.:*(α::Cyclotomic, c::T) where T<:Number = c*α
Base.:(//)(α::Cyclotomic, c::Number) = Cyclotomic(α.coeffs.//c)
Base.:(/)(α::Cyclotomic, c::Number) = Cyclotomic(α.coeffs./c)

function Base.div(α::Cyclotomic, c::Number)
    res = similar(α, first(Base.return_types(div, (eltype(α), typeof(c)))))
    return div!(res, α, c)
end

####
#   IO

function subscriptify(n::Int)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join((Char(subscript_0 + i) for i in reverse(digits(n))))
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
            coeff_str = isone(val) ? "" : "$val*"
            exp_str = isone(exp) ? "" : "^$exp"
            print(io, sign_str, coeff_str, ζ, exp_str)
        end
    end
end

end # of module Cyclotomics
