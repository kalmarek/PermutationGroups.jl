# schreier_sims implementation

function __schreier_sims_transversal(::Type{Transversal}, P::Type)
    return Transversal{AP.inttype(P),P}
end
function __schreier_sims_transversal(::Type{SchreierTransversal}, P::Type)
    return SchreierTransversal{AP.inttype(P),P,typeof(^)}
end

function schreier_sims(gens::AbstractVector{<:AbstractPermutation})
    return schreier_sims(
        __schreier_sims_transversal(Transversal, eltype(gens)),
        gens,
    )
end

function schreier_sims(
    gens::AbstractVector{<:AbstractPermutation},
    order::Integer,
)
    return schreier_sims(
        __schreier_sims_transversal(Transversal, eltype(gens)),
        gens,
        order,
    )
end

function schreier_sims(
    Tr::Type{<:AbstractTransversal{T,P}},
    gens::AbstractVector{P},
    order::Integer,
) where {T,P<:AbstractPermutation}
    @warn "schreier_sims using order is not implemented yet"
    return schreier_sims(Tr, gens)
end

function schreier_sims(
    Tr::Type{<:AbstractTransversal{T,P}},
    gens::AbstractVector{P},
) where {T,P<:AbstractPermutation}
    sc = StabilizerChain{P,Tr}()
    for s in gens
        push!(sc, s)
    end
    return sc
end

function Base.push!(stabch::StabilizerChain, g::AbstractPermutation)
    g = sift(stabch, g)
    isone(g) && return stabch

    if istrivial(stabch)
        extend_chain!(stabch, g)
    else
        extend_gens!(stabch, g)
    end
end

function sift(stabch::StabilizerChain, g::AbstractPermutation)
    if istrivial(stabch) || isone(g)
        return g
    else
        T = transversal(stabch)
        δ = point(stabch)^g
        if δ ∉ T
            return g
        else
            r = T[δ]
            if g == r # g*inv(r) == one(g)
                return one(g)
            else
                return sift(stabilizer(stabch), g * inv(r))
            end
        end
    end
end

function extend_chain!(
    stabch::StabilizerChain{P,T},
    g::AbstractPermutation,
) where {P,T}
    @assert istrivial(stabch)
    @assert !isone(g)

    # we want to modify stabch in-place, so we access the fields directly
    push!(stabch.gens, g)
    # the special transversal constructor with a single generator
    stabch.transversal = T(AP.firstmoved(g, Base.OneTo(AP.degree(g))), g, ^)
    stabch.stabilizer = typeof(stabch)() # the next stabilizer is empty

    k = length(transversal(stabch))
    if k < order(g)
        # gᵏ stabilizes point(stabch) so is a generator for stabilizer(pts)
        extend_chain!(stabch.stabilizer, g^k)
    end
    return stabch
end

function extend_gens!(
    stabch::StabilizerChain{S,T},
    g::AbstractPermutation,
) where {S,T}
    @assert !isone(g)

    # # We have two options here:
    # # a very simple version is to recompute the whole transversal,
    # # pushing the new Schreier generators to the lower level.
    # We will use it if the depth of the SchreierTransversal becomes excessive
    if T <: SchreierTransversal
        orb_dep = depth(transversal(stabch))
        orb_len = length(orbit(stabch))
        if orb_dep > 3 + ceil(sqrt(orb_len) / 2)
            push!(stabch.gens, g)
            push!(stabch.gens, prod(rand(stabch.gens, 4)))
            recompute_transversal!(stabch)

            tr = transversal(stabch)
            for s in gens(stabch)
                for δ in tr # iterating over orbit
                    schr = tr[δ] * s * inv(tr[δ^s]) # a new Schreier generator
                    isone(schr) && continue
                    push!(stabilizer(stabch), schr)
                end
            end
            return stabch
        end
    end

    # a more "advanced" version:
    # 1. move the old points in orbit with the new generator g
    # 3. add g to the generators
    # 2. push around the newly obtained points with all the generators

    tr = transversal(stabch)
    l = length(tr)

    for δ in orbit(stabch) # move the old points with the new generator
        γ = δ^g
        if γ ∉ tr
            tr[γ] = (δ, g)
        else
            s = tr[δ] * g * inv(tr[γ]) # a new Schreier generator for the stabilizer
            isone(s) && continue
            push!(stabilizer(stabch), s)
        end
    end

    push!(stabch.gens, g) # add the new generator and
    if length(tr) > l # recompute the orbit if necessary
        for (idx, δ) in enumerate(orbit(stabch))
            idx ≤ l && continue # moving only the newly added points
            for g in gens(stabch) # now with all generators
                γ = δ^g
                if γ ∉ tr
                    tr[γ] = (δ, g)
                else
                    s = tr[δ] * g * inv(tr[γ])
                    isone(s) && continue
                    push!(stabilizer(stabch), s)
                end
            end
        end
    end
    return stabch
end
