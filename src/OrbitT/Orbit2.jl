############################################################
# Linked list with last and Nothing
############################################################

mutable struct Orbit2{T, S} <: AbstractOrbit{T,S}
    pt::T # starting point
    last::Tuple{T, Tuple{S, Nothing}} # idx -> (content, nothing)
    elts::Dict{T, Tuple{S, T}} # idx -> (content, nextptr)
end

function Orbit2(pt::T, val::S) where {T,S}
    return Orbit2(pt, (pt, (val, nothing)), Dict{T, Tuple{S, T}}())
end

Orbit2(pt::T) where {T} = Orbit2(pt, nothing)

@inline Base.in(x, orb::Orbit2) = haskey(orb.elts, x) || x == orb.last[1]

@inline function Base.push!(orb::Orbit2{T,S}, tup::Tuple{T,S}) where {T,S}
    o, val = tup
    orb.elts[first(orb.last)] = (first(last(orb.last)), o)
    orb.last = (o, (val, nothing))
    return orb
end

@inline Base.push!(orb::Orbit2{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

@inline islast(orb::Orbit2, o) = o == first(orb.last)
@inline next(orb::Orbit2, s) = (last(orb.elts[s]), last(orb.elts[s]))

# o is a point in orbit, return (value, nextptr)
@inline Base.getindex(orb::Orbit2, o) = (islast(orb, o) ? first(last(orb.last)) : first(orb.elts[o]))

@inline Base.iterate(orb::Orbit2) = orb.pt, orb.pt
@inline Base.iterate(orb::Orbit2, s) = (islast(orb, s) ? nothing : next(orb, s))
@inline Base.length(orb::Orbit2) = length(orb.elts) + 1

@inline function Base.:(==)(o1::Orbit2, o2::Orbit2)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit2) = orb.pt
