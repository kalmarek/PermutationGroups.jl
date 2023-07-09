# Group Interface

Base.one(G::PermGroup{P}) where {P} = Permutation(one(P), G)
@static if VERSION < v"1.7"
    function GroupsCore.order(::Type{T}, G::AbstractPermutationGroup) where {T}
        if !isdefined(G, :order)
            ord = order(StabilizerChain(G))
            G.order = ord
        end
        return convert(T, G.order)
    end
else
    function GroupsCore.order(::Type{T}, G::AbstractPermutationGroup) where {T}
        if !isdefined(G, :order, :sequentially_consistent)
            ord = order(StabilizerChain(G))
            @atomic G.order = ord
        end
        return convert(T, G.order)
    end
end
GroupsCore.gens(G::PermGroup) = Permutation.(G.__gens_raw, Ref(G))

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{Gr},
) where {Gr<:PermGroup}
    G = rs[]
    res = Perms.perm(one(G))
    sc = StabilizerChain(G)
    for layer in sc
        pt = rand(rng, orbit(layer))
        res = transversal(layer)[pt] * res
    end
    return Permutation(res, G)
end

# GroupElement Interface
Base.parent(g::Permutation) = g.parent

function Base.:(==)(g::Permutation, h::Permutation)
    return parent(g) === parent(h) && g.perm == h.perm
end

function Base.deepcopy_internal(g::Permutation, stackdict::IdDict)
    return Permutation(Base.deepcopy_internal(g.perm, stackdict), parent(g))
end

Base.inv(g::Permutation) = Permutation(inv(g.perm), parent(g))

function Base.:(*)(g::Permutation, h::Permutation)
    parent(g) === parent(h) ||
        error("Cannot multiply elements from different permutation groups")
    return Permutation(g.perm * h.perm, parent(g))
end

### Performance optimizations

GroupsCore.order(::Type{T}, p::Permutation) where {T} = order(T, p.perm)

#### end of Group interface
