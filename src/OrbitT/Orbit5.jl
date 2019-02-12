############################################################
# linked list; last elt points to itself
############################################################

mutable struct Orbit5{T, S} <: AbstractOrbit{T,S}
    pt::T # starting point
    last::T
    elts::Dict{T, Tuple{S, T}} # idx -> (content, nextptr)
end

function Orbit5(pt::T, val::S) where {T,S}
    return Orbit5(pt, pt, Dict{T, Tuple{S, T}}(pt=>(val,pt)))
end

Orbit5(pt::T) where {T} = Orbit5(pt, nothing)

@inline Base.in(x, orb::Orbit5) = haskey(orb.elts, x)

@inline function Base.push!(orb::Orbit5{T,S}, tup::Tuple{T,S}) where {T,S}
    new, val = tup
    orb.elts[orb.last] = (orb.elts[orb.last][1], new)
    orb.last = new
    orb.elts[new] = (val, new)
    return orb
end

@inline Base.push!(orb::Orbit5{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

# o is a point in orbit, return (value, nextptr)
@inline Base.getindex(orb::Orbit5, o) = orb.elts[o][1]

@inline Base.iterate(orb::Orbit5) = orb.pt, orb.pt

@inline islast(orb::Orbit5, s) = s == orb.last
@inline next(orb::Orbit5, s) = (last(orb.elts[s]), last(orb.elts[s]))

# s is the previous elt
@inline Base.iterate(orb::Orbit5, s) = (islast(orb, s) ? nothing : next(orb, s))

@inline Base.length(orb::Orbit5) = length(orb.elts)

@inline function Base.:(==)(o1::Orbit5, o2::Orbit5)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit5) = orb.pt
