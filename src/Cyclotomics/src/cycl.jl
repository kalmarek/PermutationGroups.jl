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

# Base.@propagate_inbounds function Base.setindex!(α::Cyclotomic, val, itr)
#     for exp in itr
#         α[exp] = val
#     end
#     return itr
# end

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
