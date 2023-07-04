"""
    PermGroup(gens...)
Permutation group generated by `gens`.

`PermGroup`s are by definition sub-groups of the full symmetric group. Order
and stabilizer chain are computed (and cached) _when needed_.
"""
mutable struct PermGroup{P<:AbstractPermutation,T<:AbstractTransversal} <:
               AbstractPermutationGroup
    __gens_raw::Vector{P}
    @atomic stabchain::StabilizerChain{P,T}
    @atomic order::BigInt

    function PermGroup(
        gens::AbstractVector{<:AbstractPermutation},
        T::Type{<:AbstractTransversal},
    )
        return new{eltype(gens),T}(map(Perms.perm, gens))
    end
end

function PermGroup(gens::AbstractVector{<:AbstractPermutation})
    return PermGroup(gens, Transversal(eltype(gens)))
end

PermGroup(gens::Vararg{P,N}) where {P,N} = PermGroup(collect(gens))

struct Permutation{P,G<:PermGroup} <: AbstractPermutation
    perm::P
    parent::G
end

__gens_raw(G::PermGroup) = G.__gens_raw

Perms.perm(p::Permutation) = Perms.perm(p.perm)

# Perms.Perm interface
Perms.inttype(::Type{<:Permutation{P}}) where {P} = Perms.inttype(P)
Perms.degree(p::Permutation) = Perms.degree(p.perm)
Base.:^(n::Integer, p::Permutation) = n^p.perm
Base.one(p::Permutation) = Permutation(one(p.perm), parent(p))

# misc
"""
    StabilizerChain(G::PermGroup)
Construct the stabilizer chain for group `G`.

The first call on a particular group `G` will construct the chain from `gens(G)`
and complete it by the deterministic Schreier-Sims algorithm.
The subsequent calls just return the cached data structure.
"""
function StabilizerChain(G::PermGroup)
    if !isdefined(G, :stabchain, :sequentially_consistent)
        stabchain = schreier_sims(__gens_raw(G))
        # this may take some time, so let's check again
        if !isdefined(G, :stabchain, :sequentially_consistent)
            @atomic G.stabchain = stabchain
        end
    end
    return G.stabchain
end

basis(G::AbstractPermutationGroup) = basis(StabilizerChain(G))