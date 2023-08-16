"""
    StabilizerChain{P<:AbstractPermutation}
A structure to represent (partial) stabilizer chain.

If `G` is a group generated by a set of permutations, then for a choice of basis
> `β = (β₀, β₁, …, βₙ)`

we obtain a sequence of groups
> `G = G₀ ≥ G₁ ≥ … ≥ G₀ = {1}`.

where `Gₖ₊₁ ≤ Stab_{Gₖ}(βₖ)`. A `stbch::PointStabilizer` typically represent `Gₖ`
in the following sense:
 * `gens(stbch)` is the generating set for `Gₖ`,
 * `point(stbch)` returns `βₖ`
 * `transversal(stbch)` is the transversal of `βₖ` w.r.t. to `gens(pts)` and
 * `stabilizer(stbch)::PointStabilizer` represents `Gₖ₊₁`.

By convention `stbch` represents `G₀` (i.e. the trivial group) if
`istrivial(stbch) == true`.
"""
mutable struct StabilizerChain{P<:AbstractPermutation,T<:AbstractTransversal}
    gens::Vector{P}
    transversal::T
    stabilizer::StabilizerChain{P,T}

    StabilizerChain{P,T}() where {P,T} = new{P,T}(P[])
end

GroupsCore.gens(pts::StabilizerChain) = pts.gens
stabilizer(pts::StabilizerChain) = pts.stabilizer
transversal(pts::StabilizerChain) = pts.transversal
orbit(pts::StabilizerChain) = orbit(transversal(pts))

point(pts::StabilizerChain) = first(transversal(pts))
GroupsCore.istrivial(pts::StabilizerChain) = isempty(gens(pts))

GroupsCore.order(stabch::StabilizerChain) = order(BigInt, stabch)
function GroupsCore.order(::Type{I}, stabch::StabilizerChain) where {I}
    istrivial(stabch) && return one(I)
    return convert(I, length(transversal(stabch))) *
           order(I, stabilizer(stabch))
end

# iteration over layers of StabilizerChain
function Base.iterate(stabch::StabilizerChain, state = stabch)
    istrivial(state) && return nothing
    return state, stabilizer(state)
end

Base.length(stabch::StabilizerChain) = depth(stabch)
Base.eltype(::Type{SC}) where {SC<:StabilizerChain} = SC

__permtype(::Type{<:StabilizerChain{P}}) where {P} = P

"""
    recompute_transversal!(sc::StabilizerChain)
Recompute the Schreier tree of `sc`.

This allows shallower Schreier trees after new generators were added to `sc`.
"""
function recompute_transversal!(sc::StabilizerChain{P,T}) where {P,T}
    pt = if isdefined(sc, :transversal)
        first(transversal(sc))
    else
        firstmoved(first(gens(sc)))
    end
    sc.transversal = T(pt, gens(sc))
    return sc.transversal
end

function depth(stabch::StabilizerChain)
    depth = 0
    while !istrivial(stabch)
        depth += 1
        stabch = stabilizer(stabch)
    end
    return depth
end

basis(stabch::StabilizerChain) = [first(transversal(sc)) for sc in stabch]

function Perms.perm(
    sc::PermutationGroups.StabilizerChain{P},
    baseimages::AbstractVector{<:Integer},
) where {P}
    @assert depth(sc) == length(baseimages)
    word = Vector{P}(undef, depth(sc))

    for (idx, β, layer) in zip(eachindex(word), baseimages, sc)
        for i in Base.OneTo(idx - 1)
            l = word[i]
            β = β^l
        end
        word[idx] = inv(transversal(layer)[β])
    end
    res = inv(*(word...))
    return res
end

function sgs(stabch::StabilizerChain{S}) where {S}
    strong_gens = S[]
    while !istrivial(stabch)
        union!(strong_gens, gens(stabch))
        stabch = stabilizer(stabch)
    end
    return strong_gens
end

function __print(io::IO, sc::StabilizerChain, indent)
    println(io, indent, "┗━┳━ Stabilizer:")
    if istrivial(sc)
        print(io, indent, "  ┗━ (the trivial group)")
    else
        print(io, indent, "  ┠─┬─ ")
        println(io, "Generated by ")
        for (i, g) in enumerate(gens(sc))
            if i ≠ length(gens(sc))
                println(io, indent, "  ┃ ├─ ", g)
            else
                println(io, indent, "  ┃ └─ ", g)
            end
        end
        print(io, indent, "  ┠─── ")
        show(io, orbit(sc))
        println(io)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", sc::StabilizerChain)
    print(io, "┏ ")
    show(io, sc)
    println(io, ':')
    idnt = 0
    while true
        indent = ' '^idnt
        __print(io, sc, indent)
        istrivial(sc) && break
        sc = sc.stabilizer
        idnt += 2
    end
end

function Base.show(io::IO, sc::StabilizerChain)
    print(
        io,
        "Stabilizer chain of depth $(depth(sc)) and order $(order(sc)) with basis ",
    )
    return Base.show(io, convert(Vector{Int}, basis(sc)))
end

# iterating over the leafs of StabilizerChain tree
struct Leafs{T}
    iters::Vector{T}
    total_len::Int
end

Base.length(lfs::Leafs) = lfs.total_len
Base.eltype(::Type{<:Leafs{<:AbstractTransversal{T,S}}}) where {T,S} = S

function leafs(stabch::StabilizerChain{P,T}) where {P,T}
    transversals = T[]
    sc = stabch
    while !istrivial(sc)
        push!(transversals, transversal(sc))
        sc = stabilizer(sc)
    end

    return Leafs(transversals, prod(length, transversals; init = 1))
end

function Base.iterate(lfs::Leafs{<:AbstractTransversal})
    states = last.(iterate.(lfs.iters))

    partial_prods = map(1:length(lfs.iters)-1) do idx
        tr = lfs.iters[idx]
        return tr[first(tr)]
    end

    state = (states = states, partial_prods = partial_prods)

    res = one(eltype(lfs))
    return res, state
end

function Base.iterate(lfs::Leafs{<:AbstractTransversal}, state)
    isempty(lfs.iters) && return nothing
    next = iterate(lfs.iters[end], state.states[end])
    depth = length(lfs.iters)

    if !isnothing(next)
        pt, state.states[end] = next
        tr = lfs.iters[end]
        res = depth > 1 ? tr[pt] * state.partial_prods[end] : tr[pt]
        return res, state
    else
        while isnothing(next)
            # we need to reset the first iterator, i.e. we reached the bottom of the tree
            depth == 1 && return nothing
            # @debug "resetting $d-th iterator" state.states[d]
            _, state.states[depth] = iterate(lfs.iters[depth])
            depth -= 1
            next = iterate(lfs.iters[depth], state.states[depth])
        end
        pt, state.states[depth] = next
        tr = lfs.iters[depth]
        g = depth > 1 ? tr[pt] * state.partial_prods[depth-1] : tr[pt]

        for idx in depth:length(state.partial_prods)
            state.partial_prods[idx] = g
        end

        res = state.partial_prods[end]
        return res, state
    end
end
