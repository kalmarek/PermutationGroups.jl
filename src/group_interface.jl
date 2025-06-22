# Group Interface

Base.one(G::PermGroup{P}) where {P} = Permutation(one(P), G)

function GroupsCore.order(::Type{T}, G::AbstractPermutationGroup) where {T}
    if !isdefined(G, :order, :sequentially_consistent)
        ord = order(StabilizerChain(G))
        @atomic G.order = ord
    end
    return convert(T, G.order)
end

GroupsCore.gens(G::PermGroup) = Permutation.(G.__gens_raw, Ref(G))

function Random.Sampler(
    RNG::Type{<:Random.AbstractRNG},
    G::AbstractPermutationGroup,
    repetition::Random.Repetition = Val(Inf),
)
    return Random.SamplerTrivial(G)
end

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{Gr},
) where {Gr<:PermGroup}
    G = rs[]
    sc = StabilizerChain(G)
    P = __permtype(typeof(sc))
    word = Vector{P}(undef, length(sc))
    @inbounds for (idx, layer) in enumerate(sc)
        pt = rand(rng, transversal(layer))
        word[end-idx+1] = transversal(layer)[pt]
    end
    res = *(word...)
    return Permutation(res, G)
end

# GroupElement Interface
Base.parent(g::Permutation) = g.parent

function Base.:(==)(g::Permutation, h::Permutation)
    return parent(g) === parent(h) && g.perm == h.perm
end

function Base.deepcopy_internal(g::Permutation, id::IdDict)
    haskey(id, g) && return id[g]
    haskey(id, g.perm) && return Permutation(id[g.perm], parent(g))
    id[g] = Permutation(Base.deepcopy_internal(g.perm, id), parent(g))
    return id[g]
end

Base.inv(g::Permutation) = Permutation(inv(g.perm), parent(g))

function Base.:(*)(g::Permutation, h::Permutation)
    if parent(g) !== parent(h)
        @assert AP.perm(h) in parent(g)
    end
    return Permutation(g.perm * h.perm, parent(g))
end

### Performance optimizations

GroupsCore.order(::Type{T}, p::Permutation) where {T} = order(T, p.perm)

#### end of Group interface
