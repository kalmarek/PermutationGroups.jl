############################################################
# Linked list with last (without nothing)
############################################################

mutable struct Orbit3{T, S} <: AbstractOrbit{T,S}
    pt::T # starting point
    last::Tuple{T, S} # idx, content
    elts::Dict{T, Tuple{S, T}} # idx -> (content, nextptr)
end

function Orbit3(pt::T, val::S) where {T,S}
    return Orbit3(pt, (pt, val), Dict{T, Tuple{S, T}}())
end

Orbit3(pt::T) where {T} = Orbit3(pt, nothing)

@inline Base.in(x, orb::Orbit3) = haskey(orb.elts, x) || x == orb.last[1]

@inline function Base.push!(orb::Orbit3{T,S}, tup::Tuple{T,S}) where {T,S}
    new, val = tup
    orb.elts[first(orb.last)] = (last(orb.last), new)
    orb.last = (new, val)
    return orb
end

@inline Base.push!(orb::Orbit3{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

# o is a point in orbit, return (value, nextptr)
@inline Base.getindex(orb::Orbit3, o) = (islast(orb, o) ? last(orb.last) : first(orb.elts[o]))

@inline islast(orb::Orbit3, s) = s == first(orb.last)
@inline next(orb::Orbit3, s) = (last(orb.elts[s]), last(orb.elts[s]))

@inline Base.iterate(orb::Orbit3) = orb.pt, orb.pt
@inline Base.iterate(orb::Orbit3, s) = (islast(orb, s) ? nothing : next(orb, s))
@inline Base.length(orb::Orbit3) = length(orb.elts) + 1

@inline function Base.:(==)(o1::Orbit3, o2::Orbit3)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit3) = orb.pt
