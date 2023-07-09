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

"""
    recompute_transversal!(sc::StabilizerChain)
Recompute the Schreier tree of `sc`.

This allows shallower Schreier trees after new generators were added to `sc`.
"""
function recompute_transversal!(sc::StabilizerChain{P,T}) where {P,T}
    return sc.transversal = T(first(transversal(sc)), gens(sc))
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

struct Leafs{T}
    iters::Vector{T}
    total_len::Int
end

Base.length(lfs::Leafs) = lfs.total_len
Base.eltype(::Type{<:Leafs{<:AbstractTransversal{T,S}}}) where {T,S} = S

function leafs(stabch::StabilizerChain{P,T}) where {P,T}
    transversals = let trs = T[], sch = stabch
        while !istrivial(sch)
            push!(trs, transversal(sch))
            sch = stabilizer(sch)
        end
        trs
    end

    return Leafs(transversals, prod(length, transversals; init = 1))
end

function Base.iterate(lfs::Leafs{<:AbstractTransversal})
    states = last.(iterate.(lfs.iters))

    labels = [tr[first(tr)] for tr in Iterators.reverse(lfs.iters)]

    state = (states = states, labels = labels)
    res = *(state.labels...)
    # @info state res
    return res, state
end

function Base.iterate(lfs::Leafs{<:AbstractTransversal}, state)
    tr, st = lfs.iters[end], state.states[end]
    next = iterate(tr, st)

    if !isnothing(next)
        pt, state.states[end] = next
        state.labels[begin] = tr[pt]
        res = *(state.labels...)
        # @info state res
        return res, state
    else
        depth = length(lfs.iters)
        while isnothing(next)
            # @debug "resetting $d-th iterator" state.states[d]
            _pt, state.states[depth] = iterate(lfs.iters[depth])
            state.labels[end-depth+1] = lfs.iters[depth][_pt]
            depth -= 1
            depth == 0 && return nothing
            next = iterate(lfs.iters[depth], state.states[depth])
        end
        # @debug "advancing $d-th iterator" state.states
        # move forward the next iterator
        tr = lfs.iters[depth]
        pt, state.states[depth] = next
        state.labels[end-depth+1] = tr[pt]
        # @debug state.states

        # @show length(state.partial_products) - d
        res = *(state.labels...)
        # @info state res
        return res, state
    end
end

# function Base.iterate(stabch::StabilizerChain, state)
#     next = iterate(last(state.transversals), last(state.states))
#     @inbounds if !isnothing(next)
#         # @debug "next leaf" state.states
#         begin
#             pt, st = next
#             g = state.transversals[end][pt]
#             state.states[end] = st
#             # @debug "after" state.states
#             # we never touch state.elts[end]
#             res = state.elts[end-1] * g
#         end
#         return res, state
#     else
#         # @debug "backtrack" state.states

#         d = length(state.transversals)
#         while isnothing(next)
#             _, state.states[d] = iterate(state.transversals[d])
#             # @debug "resetting $d-th iterator" state.states[d]
#             d -= 1
#             d == 0 && return nothing
#             next = iterate(state.transversals[d], state.states[d])
#         end
#         # @debug "advancing $d-th iterator" state.states
#         # move forward the next iterator
#         pt, st = next
#         g = state.transversals[d][pt]
#         state.states[d] = st
#         # @debug state.states

#         @time state.elts[d] = d > 1 ? state.elts[d-1] * g : g
#         @show length(state.transversals) - d
#         @time begin
#             # the following elts are the same since tr[first(tr)] is always id!
#             # @views fill!(state.elts[d+1:end], state.elts[d])
#             for i in d+1:length(state.transversals)-1
#                 state.elts[i] = state.elts[i-1]
#             end
#         end

#         return state.elts[end], state
#     end
# end

# function Base.iterate(stabch::StabilizerChain{P,T}) where {P,T}
#     # initialize here
#     istrivial(stabch) && return nothing

#     transversals = let trs = T[], sch = stabch
#         while !istrivial(sch)
#             push!(trs, transversal(sch))
#             sch = stabilizer(sch)
#         end
#         trs
#     end

#     states = last.(iterate.(transversals))
#     elts = [tr[first(tr)] for tr in transversals]

#     state = (transversals = transversals, states = states, elts = elts)

#     return last(elts), state
# end
