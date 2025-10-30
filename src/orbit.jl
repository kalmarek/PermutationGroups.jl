#########################################################
# Orbit and Transversals
#########################################################
"""
    AbstractOrbit{T}
Abstract type representing abstract orbits of elements of type `T`.
"""
abstract type AbstractOrbit{T} end

struct NotInOrbit <: Exception
    x::Any
    first::Any
end
function Base.showerror(io::IO, e::NotInOrbit)
    return print(io, e.x, " was not found in the orbit of ", e.first)
end

Base.eltype(::Type{<:AbstractOrbit{T}}) where {T} = T
struct Orbit{T} <: AbstractOrbit{T}
    points::Vector{T}

    Orbit{T}(v::AbstractVector) where {T} = new{T}(v)

    function Orbit{T}(
        pt,
        gens::AbstractVector{<:GroupElement},
        op = ^,
    ) where {T}
        points = Orbit{T}(T[pt])
        for o in points
            for g in gens
                γ = op(o, g)
                if γ ∉ points
                    push!(points, γ; check = false)
                end
            end
        end
        return points
    end
end

Orbit(v::AbstractVector{T}) where {T} = Orbit{T}(v)
function Orbit(pt, gens::AbstractVector{<:GroupElement}, op = ^)
    return Orbit{typeof(pt)}(pt, gens, op)
end

Base.in(pt, orb::Orbit) = pt in orb.points
function Base.push!(orb::Orbit, pt; check = true)
    check && pt in orb && return orb
    push!(orb.points, pt)
    return orb
end

Base.length(orb::Orbit) = Base.length(orb.points)
@inline Base.iterate(orb::Orbit) = Base.iterate(orb.points)
@inline Base.iterate(orb::Orbit, state) = Base.iterate(orb.points, state)

Base.last(orb::Orbit) = last(orb.points)

function Base.show(io::IO, ::MIME"text/plain", orb::Orbit)
    return print(io, "Orbit of length $(length(orb)): ", orb.points)
end
function Base.show(io::IO, ::MIME"text/plain", orb::Orbit{<:Integer})
    return print(
        io,
        "Orbit of length $(length(orb)): ",
        convert(Vector{Int}, orb.points),
    )
end

function Base.show(io::IO, orb::Orbit)
    return print(io, typeof(orb), ':', ' ', orb.points)
end
function Base.show(io::IO, orb::Orbit{<:Integer})
    return print(io, typeof(orb), ':', ' ', convert(Vector{Int}, orb.points))
end

function Base.rand(rng::Random.AbstractRNG, st::Random.SamplerTrivial{<:Orbit})
    orb = st[]
    return rand(rng, orb.points)
end

"""
    AbstractTransversal{T,S} <: AbstractOrbit{T}
Abstract type representing the bijection of orbit oand orbit representatives.

`T` is the type of elements in the orbit, while `S` is the type of the
representatives. When `tr` is a transversal of `x` and `g` is a `GroupElement`
then `tr[x^g]` returns the representative of the `Stab(x)g` - the right coset
of the stabilizer of `x`.

## Methods to implement:
 * Constructors:
  - `Transversal(x, g::GroupElement[, action=^])` a specific constructor for a
    cyclic group
  - `Transversal(x, S::AbstractVector{<:GroupElement}[, action=^])` the default
    constructor
 * `Base.getindex(tr::T, n::Integer)` - return the coset representative
   corresponding to `n`, i.e. a group element `g` such that `first(tr)^g == n`.
   If no such element exists a `NotInOrbit` exception will be thrown.
 * Iteration protocol, iterating over points in the orbit.
"""
abstract type AbstractTransversal{T,S} <: AbstractOrbit{T} end

Base.length(tr::AbstractTransversal) = length(orbit(tr))
@inline Base.iterate(tr::AbstractTransversal) = iterate(orbit(tr))
@inline Base.iterate(tr::AbstractTransversal, state) = iterate(orbit(tr), state)
Base.last(tr::AbstractTransversal) = last(orbit(tr))

function Base.rand(
    rng::Random.AbstractRNG,
    st::Random.SamplerTrivial{<:AbstractTransversal},
)
    tr = st[]
    return rand(rng, orbit(tr))
end

function Base.getindex(tr::AbstractTransversal, pt)
    if pt in tr
        coset_representative(pt, tr)
    else
        throw(NotInOrbit(pt, first(tr)))
    end
end

struct Transversal{T,S<:GroupElement} <: AbstractTransversal{T,S}
    orbit::Orbit{T}
    representatives::Dict{T,S}

    function Transversal{T,S}(pt, g::GroupElement, op = ^) where {T,S}
        orb = Orbit{T}(T[pt])
        reps = Dict{T,S}(pt => one(g))
        tr = new{T,S}(orb, reps)

        for δ in orbit(tr)
            γ = op(δ, g)
            if γ ∉ tr
                tr[γ] = (δ, g)
            end
        end
        return tr
    end

    function Transversal{T,S}(
        pt,
        gens::AbstractVector{S},
        op = ^,
    ) where {T,S<:GroupElement}
        @assert !isempty(gens)
        orb = Orbit{T}(T[pt])
        reps = Dict{T,S}(pt => one(first(gens)))
        tr = new{T,S}(orb, reps)

        for δ in orbit(tr)
            for g in gens
                γ = op(δ, g)
                if γ ∉ tr
                    tr[γ] = (δ, g)
                end
            end
        end
        return tr
    end
end

orbit(tr::Transversal) = tr.orbit
Base.in(pt, tr::Transversal) = pt in keys(tr.representatives)
coset_representative(pt, tr::Transversal) = tr.representatives[pt]

function Base.setindex!(tr::Transversal, pt0_g::Tuple, pt1)
    @assert pt1 ∉ tr
    pt0, g = pt0_g
    push!(tr.orbit, pt1)
    tr.representatives[pt1] = tr.representatives[pt0] * g
    return tr
end

struct SchreierTransversal{T,S,Op} <: AbstractTransversal{T,S}
    orbit::Orbit{T}
    representatives::Dict{T,S}

    function SchreierTransversal{T,S,Op}(
        pt,
        g::GroupElement,
        op::Op,
    ) where {T,S,Op<:Function}
        orbit = Orbit{T}(T[pt])
        reps = Dict{T,S}(pt => one(g))
        δ = pt
        for δ in orbit
            γ = op(δ, g)
            if !(γ in keys(reps))
                push!(orbit, γ)
                reps[γ] = g
            end
        end
        return new{T,S,typeof(op)}(orbit, reps)
    end

    function SchreierTransversal{T,S,Op}(
        pt,
        gens::AbstractVector{S},
        op::Op,
    ) where {T,S<:GroupElement,Op<:Function}
        @assert !isempty(gens)
        orbit = Orbit{T}(T[pt])
        reps = Dict{T,S}(pt => one(first(gens)))
        δ = pt
        for δ in orbit
            for g in gens
                γ = op(δ, g)
                if !(γ in keys(reps))
                    push!(orbit, γ)
                    reps[γ] = g
                end
            end
        end
        return new{T,S,typeof(op)}(orbit, reps)
    end
end

function SchreierTransversal{T,S}(pt, gens, op = ^) where {T,S}
    return SchreierTransversal{T,S,typeof(op)}(pt, gens, op)
end

function SchreierTransversal{T,S,Op}(pt, gens, op = Op.instance) where {T,S,Op}
    return SchreierTransversal{T,S,typeof(op)}(pt, gens, op)
end

orbit(tr::SchreierTransversal) = tr.orbit
Base.in(pt, tr::SchreierTransversal) = pt in keys(tr.representatives)
@inline action(::SchreierTransversal{T,S,Op}) where {T,S,Op} = Op.instance

function Base.setindex!(tr::SchreierTransversal, pt0_g::Tuple, pt1)
    @assert pt1 ∉ tr
    _, g = pt0_g
    push!(tr.orbit, pt1)
    tr.representatives[pt1] = g
    return tr
end

depth(tr::SchreierTransversal) = depth(last(orbit(tr)), tr)

function depth(pt, tr::SchreierTransversal)
    depth = 1
    op = action(tr)
    while pt ≠ first(tr)
        g = tr.representatives[pt]
        pt = op(pt, inv(g))
        depth += 1
    end
    return depth
end

function coset_representative(
    pt0,
    tr::SchreierTransversal,
    depth::Integer = depth(pt0, tr);
    divide_threshold = 4,
)
    if depth ≤ divide_threshold
        # @info "base case" depth
        return coset_representative_recursive(pt0, tr)
    else
        k = depth ÷ 2
        pt1 = _descend(tr, pt0, k)
        # @info "recursive" depth k depth - k
        g0 = coset_representative(pt0, tr, k)
        g1 = coset_representative(pt1, tr, depth - k)
        return g1 * g0
    end
end

function _descend(tr::SchreierTransversal, pt, depth)
    op = action(tr)
    while pt ≠ first(tr) && depth > 1
        g = tr.representatives[pt]
        pt = op(pt, inv(g))
    end
    return pt
end

function coset_representative_recursive(
    pt0,
    tr::SchreierTransversal,
    target = first(tr),
)
    g = tr.representatives[target]
    pt0 == target && return g
    g0 = tr.representatives[pt0]
    op = action(tr)
    pt1 = op(pt0, inv(g0))
    # return coset_representative_rec(pt1^inv(g0), tr, target=target) * g0 * g
    # doing tail call optimization manually below
    if pt1 == target
        return g0 * g
    else
        g1 = tr.representatives[pt1]
        pt2 = op(pt1, inv(g1))
        if pt2 == target
            return g1 * g0 * g
        else
            g2 = tr.representatives[pt2]
            pt3 = op(pt2, inv(g2))
            if pt3 == target
                return g2 * g1 * g0 * g
            else
                return coset_representative(pt3, tr) * g2 * g1 * g0 * g
            end
        end
    end
end

## constructors for convienience
for Trans_t in (:Transversal, :SchreierTransversal)
    @eval begin
        function $Trans_t{T}(
            pt,
            gens::AbstractVector{<:GroupElement},
            op = ^,
        ) where {T}
            return $Trans_t{T,eltype(gens)}(pt, gens, op)
        end
        function $Trans_t(pt, gens::AbstractVector{<:GroupElement}, op = ^)
            return $Trans_t{typeof(pt)}(pt, gens, op)
        end
    end
end
