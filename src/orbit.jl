#########################################################
# Orbit Types and definitions
#########################################################

@inline Base.eltype(::AbstractOrbit{T}) where {T} = T

# each of them implements _Orbit Protocol_

include("OrbitT/Orbit1.jl")
include("OrbitT/Orbit2.jl")
include("OrbitT/Orbit3.jl")
include("OrbitT/Orbit4.jl")
include("OrbitT/Orbit5.jl")

const Orbit = Orbit1 # but the naive Vector+Dict is still the fastest...

for OrbT in [Symbol("Orbit", i) for i in 1:5]
    @eval begin
        @doc doc"""
            Orbit(gens::Vector{<:GrouElem}, pt[, ^])
        Compute the orbit of a point `pt` under the action of group generated by `gens`.

        Only the consecutive points are stored in the orbit.
        """
        function $OrbT(gens::AbstractVector{<:GroupElement}, pt::T, op = ^) where {T}
            orb = $OrbT(pt)
            for o in orb
                for g in gens
                    γ = op(o, g)
                    if γ ∉ orb
                        push!(orb, γ)
                    end
                end
            end
            return orb
        end

        @doc doc"""
            Orbit(::Type{Transversal}, gens::Vector{<:GrouElem}, pt[, ^])
        Compute the orbit of a point `pt` under the action of group generated by `gens`.

        Along the points the right transversal is stored as explicit group elements.
        """
        function $OrbT(
            ::Type{Transversal},
            gens::Vector{<:GroupElement},
            pt::T,
            op = ^,
        ) where {T}
            orb = $OrbT(pt, one(parent(first(gens))))
            for o in orb
                for g in gens
                    γ = op(o, g)
                    if γ ∉ orb
                        push!(orb, (γ, orb[o] * g))
                    end
                end
            end
            return orb
        end
    end
end
