(G::Group)(gens::AbstractVector{<:Perm}) = PermGroup(gens)

# embeddings

function emb!(res::AP, g::AP, V::AbstractVector{<:Integer}) where {AP<:AbstractPerm}
    @assert degree(res) >= degree(g)
    res.d .= Base.OneTo(degree(res))
    res.d[V] .= g.d
    return res
end

emb(g::Perm, n::Integer) = emb!(Perm(n), g, 1:degree(g))

function (G::PermGroup)(p::AbstractPerm)
    g = degree(p) == degree(G) ? p : emb(p, degree(G))
    return Permutation(g, G)
end

@doc doc"""
    base(G::PermGroup)
Get the base of the group `G`.
"""
base(G::AbstractPermutationGroup) = base(StabilizerChain(G))

@doc doc"""
    sgs(G::PermGroup)
Compute strong generating set for group `G`.

If a stabilizer chain for `G` has been already stored in `G`, `sgs` will be
just retrieved from it.
"""
sgs(G::AbstractPermutationGroup) = sgs(StabilizerChain(G))

@doc doc"""
    in(g::perm, G::PermGroup)
Membership test for permutation group `G` by `sift`ing `g` through `StabilizerChain(G)`.
"""
function Base.in(g::Perm, G::PermGroup)
    g = degree(g) < degree(G) ? emb(g, degree(G)) : g

    sc = StabilizerChain(G)
    h, depth = sift(g, sc)
    depth ≤ length(sc) && return false
    return ifelse(isone(h), true, false)
end

@doc doc"""
	transversals(G::PermGroup)
Return the transversals (as a Vector) of a permutation group `G`.
"""
transversals(G::AbstractPermutationGroup) = transversals(StabilizerChain(G))

@doc doc"""
    perm_by_baseimages(G::PermGroup, base_images::AbstractVector{<:Integer})
Return the unique permutation in `G` determined by to `base_images`
(with respect to `base(G)`).
"""
function perm_by_baseimages(
    G::AbstractPermutationGroup,
    baseimages::AbstractVector{<:Integer},
)
    @boundscheck length(baseimages) == length(base(G))
    trans = transversals(G)
    res = one(G)

    for (schr, bi) in zip(trans, baseimages)
        gr_word = Word(schr.gens_inv, schr.orb, bi, schr.op)
        res = gr_word(schr.gens_inv, res)
        # res = mul!(res, res, t[bi])
    end

    return G(inv(res))
end
