"""
    parity(g::AbstractPermutation)

Return the parity of number of factors in factorization of `g` into transpositions.

Return `1`` if the number is odd and `0` otherwise.
"""
parity(σ::AbstractPermutation) = __parity_generic(σ)
parity(cd::CycleDecomposition) = Int(isodd(count(iseven ∘ length, cd)))

function __parity_generic(σ::AbstractPermutation)
    to_visit = trues(degree(σ))
    parity = false
    k = 1
    @inbounds while any(to_visit)
        k = findnext(to_visit, k)
        to_visit[k] = false
        next = k^σ
        while next != k
            parity = !parity
            to_visit[next] = false
            next = next^σ
        end
    end
    return Int(parity)
end

"""
    sign(g::AbstractPermutation)

Return the sign of a permutation.

`sign` represents the homomorphism from the permutation group to the unit group
of `ℤ` whose kernel is the alternating group.
"""
Base.sign(σ::AbstractPermutation) = (-1)^parity(σ)

"""
    permtype(g::AbstractPermutation)

Return the type of permutation `g`, i.e. lengths of disjoint cycles in cycle
decomposition of `g`.

The lengths are sorted in decreasing order by default. `permtype(g)` fully
determines the conjugacy class of `g`.
"""
function permtype(σ::AbstractPermutation)
    return sort!([length(c) for c in cycles(σ) if length(c) > 1]; rev = true)
end

function GroupsCore.order(::Type{T}, σ::AbstractPermutation) where {T}
    isone(σ) && return one(T)
    return GroupsCore.order(T, cycles(σ))
end

function GroupsCore.order(::Type{T}, cd::CycleDecomposition) where {T}
    return convert(T, mapreduce(length, lcm, cd; init = 1))
end

function firstmoved(σ::AbstractPermutation)
    @inbounds for i in Base.OneTo(degree(σ))
        if i^σ ≠ i
            return i
        end
    end
    return nothing
end

function fixedpoints(p::AbstractPermutation, range = Base.OneTo(degree(p)))
    return [i for i in range if i^p == i]
end
function nfixedpoints(p::AbstractPermutation, range = Base.OneTo(degree(p)))
    return count(i -> i^p == i, range)
end

# function fixes(p::AbstractPermutation, v::AbstractVector, op = ^)
#     return all(i -> v[i] == v[op(i, p)], Base.OneTo(degree(p)))
# end

#=
function power_by_cycles(σ::AbstractPermutation, n::Integer)
    if n < 0
        return inv(σ)^-n
    elseif n == 0
        return one(σ)
    elseif n == 1
        return copy(σ)
    elseif n == 2
        return σ * σ
    elseif n == 3
        return σ * σ * σ
    elseif n == 4
        return σ * σ * σ * σ
    else
        img = Vector{inttype(σ)}(undef, deg)
        @inbounds for cycle in cycles(σ)
            l = length(cycle)
            k = n % l
            for (idx, j) in enumerate(cycle)
                idx += k
                idx = (idx > l ? idx - l : idx)
                img[j] = cycle[idx]
            end
        end
        return typeof(P)(img, false)
    end
end
=#
