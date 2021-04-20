(G::Group)(gens::AbstractVector{<:Perm}) = PermGroup(gens)

# embeddings

function emb!(res::Perm, g::Perm, V::AbstractVector{<:Integer})
    @assert degree(res) >= degree(g)
	for (idx, v) in pairs(V)
		res[v] = idx^g
	end
    return res
end

emb(g::Perm, n::Integer) = emb!(Perm(n), g, eachindex(g))

function (G::PermGroup)(p::AbstractPerm)
    g = degree(p) == degree(G) ? p : emb(p, degree(G))
    return Permutation(perm(g), G)
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
function Base.in(g::AbstractPerm, G::AbstractPermutationGroup)
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
	inverse = true,
	tmp_word = Word(eltype(eltype(G))[])
)
    @boundscheck length(baseimages) == length(base(G))
    trans = transversals(G)
    res = one(G)

	for (i, schr) in enumerate(trans)
		bi = baseimages[i]
	    bi == first(schr) && continue
	    word = append!(tmp_word, schr.gens_inv, schr.orb, bi, schr.op)
	    res = tmp_word(schr.gens_inv, res)
		empty!(tmp_word)
	end

    return inverse ? inv(res) : res
end

###

struct BaseImagesIter{I, T}
	images::Vector{I}
	transversals::T
end

Base.eltype(::Type{<:BaseImagesIter{I}}) where {I} = Vector{I}
Base.length(bitr::BaseImagesIter) = prod(length, bitr.transversals)

base_images(G::PermGroup) = BaseImagesIter(copy(base(G)), transversals(G))

@inline Base.iterate(bitr::BaseImagesIter) = bitr.images, (count=1, total=length(bitr))

@inline function Base.iterate(bitr::BaseImagesIter, state)
	state.count == state.total && return nothing
	img = next!(bitr.images, bitr.transversals)
	return img, (count=state.count+1, total=state.total)
end

@inline function next!(baseimages::AbstractVector{<:Integer}, transversals)
	for position in length(baseimages):-1:1
		t = transversals[position]
		if islast(t, baseimages[position])
			# @debug "last point in orbit: spilling at" position collect(t), baseimages[position]
			baseimages[position] = first(t)
		else
			# @debug "next point in orbit: incrementing at" position collect(t), baseimages[position]
			baseimages[position] = _next(t, baseimages[position])
			break
		end
	end
	return baseimages
end

function _next(itr, elt)
	x, state = iterate(itr)
	while x != elt
		x, state = iterate(itr, state)
	end
	return first(iterate(itr, state))
end
