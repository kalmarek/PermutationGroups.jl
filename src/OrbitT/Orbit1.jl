############################################################
# Naive Vector&Dict implementation (as fast as manual loop),
############################################################
struct Orbit1{T, S} <: AbstractOrbit{T,S}
    elts::Vector{T}
    vals::Dict{T, S}
end
Orbit1(pt::T, val::S) where {T,S} = Orbit1([pt], Dict(pt=>val))
Orbit1(pt::T) where T = Orbit1(pt, nothing)
# Orbit1(::Type{T}, ::Type{S}) where {T,S} = Orbit1(T[], Dict{T,S}())

function Orbit1(v::AbstractVector{T}) where T
    return Orbit1(Vector(v), Dict(i=>nothing for i in v))
end

@inline Base.in(x, orb::Orbit1) = haskey(orb.vals, x)
@inline function Base.push!(orb::Orbit1{T,S}, tup::Tuple{T,S}) where {T,S}
    push!(orb.elts, tup[1])
    orb.vals[tup[1]] = tup[2]
    return orb
end

@inline function Base.push!(orb::Orbit1{T,Nothing}, x::T) where {T}
    push!(orb.elts, x)
    orb.vals[x] = nothing
    return orb
end

@inline Base.getindex(orb::Orbit1, o) = orb.vals[o]

@inline Base.iterate(orb::Orbit1) = iterate(orb.elts)
@inline Base.iterate(orb::Orbit1, s) = iterate(orb.elts, s)
@inline Base.length(orb::Orbit1) = length(orb.elts)
@inline Base.:(==)(o1::Orbit1, o2::Orbit1) =
    o1.elts == o2.elts && o1.vals == o2.vals

@inline Base.first(orb::Orbit1) = orb.elts[1]
