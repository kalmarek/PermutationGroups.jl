
@doc """doc
    AbstractOrbit{T,S}
> _Orbit protocol_ consists of
> 1. constructors:
  * `OrbT(pt::T)` → `OrbT{T, Nothing}`
  * `OrbT(pt::T, val::S)` -> `OrbT{T, S}`
> 2. `Base.in(x::T, orb::OrbT{T})` → checks if `x` is in the orbit
> 3. `Base.push!` adds points to orbits:
  * `push!(::OrbT{T, Nothing}, pt) → add just a `pt` to plain orbit
  * `push!(::OrbT, (pt, val))` → add both `pt` and the corresponding value `val`)
> 4. `Base.getindex(::OrbT{T,S}, pt::T)::S` → returns the value added along `pt`
> 5. `Base.iterate` → iterates over elements in orbit
> 6. `Base.length` → the length of the orbit
> 7. `Base.eltype(::OrbT{T})` → returns element type of orbit (implemented for AbstractOrbit{T}
> 8. `Base.first` → return the initial point.
"""
abstract type AbstractOrbit{T,S} end
struct Transversal end

struct Word{I}
    gens_ptrs::Vector{I}
end

# Almost implements the OrbitProtocol.

@doc doc"""
    Schreier(gens_inv, orb[, op=^]) <: AbstractOrbit{<:GroupElem, <:Integer}
> Data structure which represents the Schreier tree, i.e. stores the orbit in
> a tree, where each branch and labeled by the inverse of a generator (i.e. is
> oriented towards the the base point).
> Schreier tree stores the inverses of generators and allows indexing by orbit
> elements to query for the stabilizer cosets representatives. That is if
> `s = Schreier(gens, pt)`, then `uᵦ = s[i]` is a group element which takes
> `pt` to `i`, i.e. `pt^uᵦ = i`. Note that if the inverse of `uᵦ` is needed
> a faster metod `getinv(schr, i)` is also available.
>
> `Schreier` almost implements the Orbit Protocol, i.e. it extends the
> following methods:
> `in`, `iterate`, `length`, `eltype`, `first`. By its definition only
> `push!(::Schreier, ::Tuple{S,I}) where I<:Integer` is supported.
"""
struct Schreier{GEl<:GroupElem, S, I<:Integer, Orb<:AbstractOrbit{S, I}, Op} <: AbstractOrbit{GEl, I}
    gens_inv::Vector{GEl}
    orb::Orb
    op::Op

    function Schreier(gens_inv::Vector{GEl}, orb::Orb, op::Op=^) where {S, GEl<:GroupElem, I<:Integer, Orb<:AbstractOrbit{S, I}, Op}
        return new{GEl, S, I, Orb, Op}(gens_inv, orb, op)
    end
end

@doc doc"""
    StabilizerChain(base, sgs, transversals)
> `StabilizerChain` struct with fields
> * `base::Vector{<:Integer}` →  stores the base of the chain
> * `sgs::Vector{Vector{<:GroupElem}}` → for each of base element the strong generating set of its stabilizer
> * `Transversals::Vector{<:Schreier}` → for each of base element the Schreier Tree of its orbit
"""
struct StabilizerChain{I, GEl, Schr <: Schreier}
    base::Vector{I}
    sgs::Vector{Vector{GEl}}
    transversals::Vector{Schr}
end

@doc doc"""
    PermGroup(gens[, stabchain])
> Permutation group (i.e. a sub-group of the full symmetric group).
> If stabilizer chain is not provided, then it will be recomputed _when needed_.
"""
mutable struct PermGroup{I<:Integer, SC<:StabilizerChain} <: AbstractAlgebra.AbstractPermutationGroup
    gens::Vector{Generic.Perm{I}}
    stabchain::SC

    function PermGroup(gens::Vector{Generic.Perm{I}}) where I<:Integer
        maxdegree = maximum(degree.(gens))
        new_gens = Generic.emb.(gens, maxdegree)
        sc = Schreier([first(new_gens)], I(1), ^)
        return new{I,
            StabilizerChain{I, Perm{I}, typeof(sc)}}(new_gens)
    end

    function PermGroup(gens::Vector{Generic.Perm{I}}, sc::SC) where {I<:Integer, SC<:StabilizerChain}
        return new{I, SC}(gens, sc)
    end
end

abstract type AbstractClassFunction{T} end # <: AbstractVector{T} ??

mutable struct Character{T, CC} <: AbstractClassFunction{T}
    vals::Vector{T}
    cc::Vector{CC}
    inv_of::Vector{Int16}

    function Character(v::AbstractVector{T}, ccls::AbstractVector{CCl},
        inv_of=_inv_of(Int16, ccls)) where {T, CCl}

        χ = new{T, CCl}(v, ccls, inv_of)

        # initial normalization
        R = parent(first(v))
        id = one(first(first(ccls)))
        if !isone(χ(id))
            χ.vals .*= inv(χ(id))
        end

        # computing the degree of χ:
        deg_χ = sqrt(inv(dot(χ, χ)))
        @debug χ.vals, deg_χ

        # renormalizing χ
        χ.vals .*= deg_χ

        return χ
    end
end

struct CCMatrix{T, C} <: AbstractMatrix{T} # M_r
    cc::Vector{C} # vector of conjugacy classes to fix the order
    r::Int # the index of conjugacy class
    m::Matrix{T} # cache of class coefficients

    function CCMatrix(cc::A, r::Int, T::Type=Int) where {C, A<:AbstractVector{C}}
        M = -ones(T, length(cc), length(cc))
        new{T, C}(cc, r, M)
    end
end

struct EigenSpaceDecomposition{R <: RingElement, M <: MatrixElem}
    basis::M
    eigspace_ptrs::Vector{Int}

    function EigenSpaceDecomposition(
        basis::MatrixElem{R},
        eigspace_ptrs::AbstractVector{<:Integer}) where R <: RingElement

        @assert eigspace_ptrs[1] == 1
        @assert eigspace_ptrs[end] == _dim(basis)+1

        return new{R, typeof(basis)}(basis, eigspace_ptrs)
    end
end
