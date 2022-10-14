@doc doc"""
    schreier_sims!(sc::StabilizerChain)
Complete the `sc`, i.e. run the _full, deterministic_ Schreier-Sims algorithm on `sc`.
"""
function schreier_sims!(sc::StabilizerChain, _order::T=0) where T<:Integer
    depth = length(sc)

    order_is_known = _order > 0

    while depth ≥ 1
        @label outer_loop

        _, S, Δ = sc[depth]
        @debug "going through orbit $(collect(Δ))"
        for β in Δ
            uᵦ = Δ[β]
            @debug "Next point in orbit" β uᵦ
            for g in S
                @debug "Scheier generator for" g
                schreier_generator = uᵦ*g*getinv(Δ, β^g)
                isone(schreier_generator) && continue
                found_new_generator = false
                @debug "a new non-trivial schreier generator" schreier_generator
                h, new_depth = sift(schreier_generator, sc)
                @debug "sifts to $h at depth $new_depth"
                if new_depth ≤ length(sc)
                    @debug "its a new strong generator at depth $new_depth !"
                    found_new_generator = true
                elseif !isone(h)
                    @debug "$new_depth > $(length(sc)): extending the chain..."
                    found_new_generator = true
                    push!(sc, firstmoved(h))
                    @debug "New base" sc.base
                end

                if found_new_generator
                    for l in (depth+1):new_depth # h fixes all points <= depth
                        @debug "recomputing Schreier trees at level $l ∈ $((depth+1):new_depth)"
                        push!(sc, h, l, recompute=true) # recompute Schreier trees
                        order
                        @debug "...done! new transversal is" sc.transversals[l]
                    end

                    if order_is_known
                        if order(T, sc) == _order
                            @goto chain_completed
                        end
                    end
                    @debug "restarting the procedure at depth" new_depth
                    depth = new_depth
                    @goto outer_loop
                end
            end
        end
        depth -= 1
        @debug "Checked all gens, decreasing depth to" depth
    end

    @label chain_completed
    order_is_known && @assert order(sc) == _order
    return sc
end

@doc doc"""
    schreier_sims(gens::Vector{perm}, B::AbstractVector{<:Integer}) → StabilizerChain
Complete `gens` to strong generator set by including `B` as (partial) base.
"""
function schreier_sims(
    gens::AbstractVector{<:AbstractPerm},
    B::AbstractVector{<:Integer}=eltype(eltype(gens))[],
    _order=0)
    sc = StabilizerChain(gens, B)
    schreier_sims!(sc, _order)
    return sc.base, sc.sgs, sc.transversals
end
