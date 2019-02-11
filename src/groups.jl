function Base.show(io::IO, PG::PrmGroup)
    init = isdefined(PG, :stabchain) ? " with stabiliser chain of length $(length(PG.stabchain))" : ""
    println(io, "Permutation group$init")
    println(io, "⟨"*join(gens(G), ", ")*"⟩")
end

@doc doc"""
    gens(G::PrmGroup)
> The original generators of the permutation group `G`.
"""
AbstractAlgebra.gens(G::PrmGroup) = G.gens

@doc doc"""
    StabilizerChain(G::PrmGroup)
> Construct the stabilizer chain for group `G`. If constructed for the first time
> the chain will be completed by the full deterministic Schreier-Sims algorithm.
> The `StabilizerChain` struct will be stored in `G` for future uses.
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
    h, depth = sift(g, sc)
    depth ≤ length(sc) && return false
    return ifelse(isone(h), true, false)
end
