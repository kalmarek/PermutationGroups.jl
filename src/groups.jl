function Base.show(io::IO, G::PrmGroup)
    init = isdefined(G, :stabchain) ? " of order $(order(StabilizerChain(G)))" : ""
    println(io, "Permutation group on $(length(gens(G))) generators$init")
    println(io, "⟨"* gensstring(gens(G), width=240)*"⟩")
end

@doc doc"""
    gens(G::PrmGroup)
> The original generators of the permutation group `G`.
"""
AbstractAlgebra.gens(G::PrmGroup) = G.gens

@doc doc"""
    StabilizerChain(G::PrmGroup)
> Construct the stabilizer chain for group `G`. The first call on a particular
> group `G` will construct the chain from `gens(G)` and complete by the
> deterministic Schreier-Sims algorithm. The subsequent calls just return the
> cached copy.
"""
function StabilizerChain(G::PrmGroup)
    if !(isdefined(G, :stabchain))
        G.stabchain = StabilizerChain(gens(G))
        schreier_sims!(G.stabchain)
    end
    return G.stabchain
end

@doc doc"""
    base(G::PrmGroup)
> Get the base of the group `G`.
"""
base(G::PrmGroup) = StabilizerChain(G).base

@doc doc"""
    sgs(G::PrmGroup)
> Compute strong generating set for group `G`. If a stabilizer chain for `G`
> has been already stored in `G`, `sgs` will be just retrieved from it.
"""
sgs(G::PrmGroup) = sgs(StabilizerChain(G))

@doc doc"""
    order(G::PrmGroup)::BigInt
> Compute the order of permutation group `G` by computing stabilizer chain. If a
> stabiliser chain for `G` has been already stored in `G`, it will be used for the
> computation.
"""
AbstractAlgebra.order(G::PrmGroup) = order(StabilizerChain(G))

@doc doc"""
    in(g::perm, G::PrmGroup)
> Membership test for permutation group `G` by `sift`ing `g` through `StabilizerChain(G)`.
"""
function Base.in(g, G::PrmGroup)
    sc = StabilizerChain(G)
    n = degree(first(gens(G)))
    if degree(g) < n
        g = Generic.emb(g, n)
    end
    h, depth = sift(g, sc)
    depth ≤ length(sc) && return false
    return ifelse(isone(h), true, false)
end
