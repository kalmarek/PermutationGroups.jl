# schreier_sims implementation

function schreier_sims(gens::AbstractVector{<:AbstractPermutation})
    return schreier_sims(Transversal(eltype(gens)), gens)
end

function schreier_sims(
    gens::AbstractVector{<:AbstractPermutation},
    order::Integer,
)
    return schreier_sims(Transversal(eltype(gens)), gens, order)
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
        x = point(stabch)
        δ = x^g
        if δ ∉ T
            return g
        else
            r = T[δ]
            g = g * inv(r)
            @assert x^g == x
            return isone(g) ? g : sift(stabilizer(stabch), g)
        end
    end
end

function extend_chain!(
    stabch::StabilizerChain{P,T},
    g::AbstractPermutation,
) where {P,T}
    @assert !isone(g)

    # we want to modify stabch in-place, so we access the fields directly
    push!(stabch.gens, g)
    # the special transversal constructor with a single generator
    stabch.transversal = T(firstmoved(g), g, ^)
    stabch.stabilizer = StabilizerChain{P,T}() # the next stabilizer is empty

    k = length(transversal(stabch))
    if k < order(g)
        # gᵏ stabilizes point(stabch) so is a generator for stabilizer(pts)
        extend_chain!(stabch, g^k)
    end
    return stabch
end

function extend_gens!(
    stabch::StabilizerChain{S,T},
    g::AbstractPermutation,
) where {S,T}
    @assert !isone(g)

    # a very simple version
    push!(stabch.gens, g)
    stabch.transversal = T(point(stabch), gens(stabch), ^)

    tr = transversal(stabch)
    for s in gens(stabch)
        for δ in tr # iterating over orbit
            schr = tr[δ] * s * inv(tr[δ^s]) # a new Schreier generator
            isone(schr) && continue
            @assert point(stabch)^schr == point(stabch)
            push!(stabilizer(stabch), schr)
        end
    end
    return stabch
end
