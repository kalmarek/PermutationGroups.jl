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
    orb.elts[orb.last[1]] = (orb.last[2], new)
    orb.last = (new, val)
    return orb
end

@inline Base.push!(orb::Orbit3{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

# o is a point in orbit, return (value, nextptr)
@inline Base.getindex(orb::Orbit3, o) = (o == orb.last[1] ? orb.last[2] : orb.elts[o][1])

@inline Base.iterate(orb::Orbit3) = orb.pt, orb.pt
# s is the previous elt
@inline Base.iterate(orb::Orbit3, s) =
    (s == orb.last[1] ? nothing : (orb.elts[s][2], orb.elts[s][2]))
@inline Base.iterate(orb::Orbit3, ::Nothing) = nothing
@inline Base.length(orb::Orbit3) = length(orb.elts) + 1

@inline function Base.:(==)(o1::Orbit3, o2::Orbit3)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit3) = orb.pt
