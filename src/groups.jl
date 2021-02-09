function Base.show(io::IO, G::PermGroup)
    init = isdefined(G, :stabchain) ? " of order $(order(StabilizerChain(G)))" : ""
    println(io, "Permutation group on $(length(gens(G))) generators$init")
    println(io, "⟨"* gensstring(gens(G), width=240)*"⟩")
end

@doc doc"""
    gens(G::PermGroup)
The original generators of the permutation group `G`.
"""
AbstractAlgebra.gens(G::PermGroup) = G.gens

@doc doc"""
    StabilizerChain(G::PermGroup)
Construct the stabilizer chain for group `G`.

The first call on a particular group `G` will construct the chain from `gens(G)`
and complete by the deterministic Schreier-Sims algorithm.
The subsequent calls just return the cached copy.

!!! Warning !!! It is users responsibility to ensure that the cached copy is
completed by call to `schreier_sims!` if the `base` or `gens` are changed
manually after groups creation.
"""
function StabilizerChain(G::PermGroup)
    if !(isdefined(G, :stabchain))
        G.stabchain = schreier_sims!(StabilizerChain(gens(G)))
    end
    return G.stabchain
end

@doc doc"""
    base(G::PermGroup)
Get the base of the group `G`.
"""
base(G::PermGroup) = base(StabilizerChain(G))

@doc doc"""
    sgs(G::PermGroup)
Compute strong generating set for group `G`.

If a stabilizer chain for `G` has been already stored in `G`, `sgs` will be
just retrieved from it.
"""
sgs(G::PermGroup) = sgs(StabilizerChain(G))

@doc doc"""
    order(G::PermGroup)::BigInt
Compute the order of permutation group `G` by computing stabilizer chain.

If a stabiliser chain for `G` has been already stored in `G`, it will be used
for the computation.
"""
AbstractAlgebra.order(G::PermGroup) = order(BigInt, StabilizerChain(G))
AbstractAlgebra.order(::Type{T}, G::PermGroup) where T =
	order(T, StabilizerChain(G))

@doc doc"""
    in(g::perm, G::PermGroup)
Membership test for permutation group `G` by `sift`ing `g` through `StabilizerChain(G)`.
"""
function Base.in(g, G::PermGroup)
    sc = StabilizerChain(G)
    n = degree(first(gens(G)))
    if degree(g) < n
        g = Generic.emb(g, n)
    end
    h, depth = sift(g, sc)
    depth ≤ length(sc) && return false
    return ifelse(isone(h), true, false)
end

@doc doc"""
	transversals(G::PermGroup)
Return the transversals (as a Vector) of a permutation group `G`.
"""
transversals(G::PermGroup) = transversals(StabilizerChain(G))

@doc doc"""
    perm_by_baseimages(G::PermGroup, base_images::AbstractVector{<:Integer})
Return the unique permutation in `G` determined by to `base_images`
(with respect to `base(G)`).
"""
function perm_by_baseimages(G::PermGroup, baseimages::AbstractVector{<:Integer})
    @boundscheck length(baseimages) == length(base(G))
	trans = transversals(G)
	res = one(G)

	for (schr, bi) in zip(trans, baseimages)
		gr_word = Word(schr.gens_inv, schr.orb, bi, schr.op)
	    res = gr_word(schr.gens_inv, res)
		# res = mul!(res, res, t[bi])
	end

    return inv(res)
end

@doc doc"""
	degree(G::PermGroup{I})::I where I
Return the degree of `G`, i.e. the length of the storage of permutations in `G`.
!!! This is an implementation detail and due to change without notice!!!
"""
Generic.degree(G::PermGroup{I}) where I = I(degree(first(gens(G))))

@doc doc"""
	one(G::PermGroup)
Return the identity of a permutation group `G`.
"""
Base.one(G::PermGroup) = Perm(degree(G))

####################################
### iteration protocol for PermGroups

function Base.iterate(G::PermGroup)
	return one(G), (deepcopy(base(G)), 1, order(G))
end

function Base.iterate(G::PermGroup, state)
	base_im, count, ord_G = state
	count == ord_G && return nothing

	basis_im = next!(base_im, transversals(G))
	g = perm_by_baseimages(G, base_im)

	return (g, (base_im, count+1, ord_G))
end

Base.eltype(::Type{<:PermGroup{I}}) where I = Perm{I}
Base.length(G::PermGroup) = order(G)
Base.size(G::PermGroup) = (length(G),)
