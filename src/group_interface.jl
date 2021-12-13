# Group Interface

Base.one(G::PermGroup) = Permutation(Perm(degree(G)), G)
GroupsCore.order(::Type{T}, G::AbstractPermutationGroup) where {T} =
    order(T, StabilizerChain(G))
GroupsCore.gens(G::PermGroup) = Permutation.(gens_raw(G), Ref(G))

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{Gr}
) where {Gr <: PermGroup}
    G = rs[]
    tr = transversals(G)
    # Ref(rng) is needed for julia-1.3
    img =rand.(Ref(rng), tr)
    return perm_by_baseimages(G, img)
end
### iteration protocol for PermGroups

Base.eltype(::Type{GT}) where {I,GT<:PermGroup{I}} = Permutation{I,GT}
Base.IteratorSize(::Type{<:AbstractPermutationGroup}) = Base.HasLength()

struct PermGroupIter{V, I, S, W}
    base_images::V
    itr::I
    state::S
    tmp_word::W
end

function Base.iterate(G::PermGroup)
    itr = base_images(G)
    tmp_word = Word(eltype(eltype(G))[])
    img, st = iterate(itr)
    permiter = PermGroupIter(img, itr, st, tmp_word)
    g = perm_by_baseimages(G, img, false, permiter.tmp_word)
    return g, (permiter, 1)
end

function Base.iterate(G::PermGroup, state)
    permiter, count = state
    count >= length(G) && return nothing
    img, st = iterate(permiter.itr, permiter.state)
    g = perm_by_baseimages(G, img, false, permiter.tmp_word)
    return g, (permiter, count+1)
end

# GroupElement Interface
Base.parent(g::Permutation) = g.parent

Base.:(==)(g::Permutation, h::Permutation) = parent(g) === parent(h) && g.perm == h.perm

Base.deepcopy_internal(g::Permutation, stackdict::IdDict) =
    Permutation(Base.deepcopy_internal(g.perm, stackdict), parent(g))

Base.inv(g::Permutation) = Permutation(inv(g.perm), parent(g))

function Base.:(*)(g::Permutation, h::Permutation)
    parent(g) === parent(h) ||
        error("Cannot multiply elements from different permutation groups")
    return Permutation(g.perm * h.perm, parent(g))
end

### Performance optimizations

Base.similar(p::Permutation{T}, ::Type{S} = T) where {T,S} =
    Permutation(similar(perm(p), T), parent(p))
Base.isone(p::AbstractPerm) = nfixedpoints(p) == degree(p)
Base.:(^)(p::Permutation, n::Integer) = parent(p)(p.perm^n)
GroupsCore.order(::Type{T}, p::Permutation) where {T} = order(T, p.perm)
Base.hash(p::Permutation, h::UInt) = hash(typeof(p), hash(p.perm, hash(parent(p), h)))

# Mutable API

function GroupsCore.one!(g::AbstractPerm)
    for i in eachindex(g)
        g[i] = i
    end
    return g
end

Base.@propagate_inbounds function GroupsCore.mul!(
    out::AbstractPerm,
    g::AbstractPerm,
    h::AbstractPerm,
)
    out = (out === h ? similar(out) : out)
    @boundscheck @assert degree(out) >= max(degree(g), degree(h))

    @inbounds for i in eachindex(out)
        out[i] = (i^g)^h #h[g[i]]
    end

    return out
end

#### end of Group interface

# accessors
gens_raw(G::PermGroup) = G.gens

# AbstractPerm Interface??
"""
	perm(g::AbstractPerm)::Perm
Return an instance of `Perm`, a parent-less permutation backed by simple vector storage.
"""
perm(g::Permutation) = g.perm
perm(g::AbstractPerm) = Perm([i^g for i in eachindex(g)], false)

Base.setindex!(g::Permutation, v::Integer, n::Integer) = g.perm[n] = v

@doc doc"""
	degree(G::PermGroup)
Return the degree of `G`, i.e. the length of the storage of permutations in `G`.
"""
degree(G::PermGroup) = return G.deg
degree(G::SymmetricGroup) = G.n
degree(p::Permutation) = degree(parent(p))
degree(p::Perm{I}) where {I} = I(length(p.d))

permtype(p::Permutation) = permtype(perm(p))
Base.sign(p::Permutation) = sign(perm(p))
