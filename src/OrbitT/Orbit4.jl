############################################################
# Linked list with Tuple{S, Union{T, Nothing}}
############################################################

mutable struct Orbit4{T, S} <: AbstractOrbit{T,S}
    pt::T # starting point
    last::T
    elts::Dict{T, Tuple{S, Union{T, Nothing}}} # idx -> (content, nextptr)
end

function Orbit4(pt::T, val::S) where {T,S}
    return Orbit4(pt, pt, Dict{T, Tuple{S, Union{T,Nothing}}}(pt=>(val,nothing)))
end

Orbit4(pt::T) where {T} = Orbit4(pt, nothing)

@inline Base.in(x, orb::Orbit4) = haskey(orb.elts, x)

@inline function Base.push!(orb::Orbit4{T,S}, tup::Tuple{T,S}) where {T,S}
    new, val = tup
    orb.elts[orb.last] = (first(orb.elts[orb.last]), new)
    orb.last = new
    push!(orb.elts, new => (val, nothing))
    return orb
end

@inline Base.push!(orb::Orbit4{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

@inline Base.getindex(orb::Orbit4, o) = first(orb.elts[o])

@inline islast(orb::Orbit4, s) = s == orb.last
@inline function next(orb::Orbit4{T}, s) where {T}
    return (last(orb.elts[s]), last(orb.elts[s]))::Tuple{T,T}
end

@inline Base.iterate(orb::Orbit4) = orb.pt, orb.pt
@inline Base.iterate(orb::Orbit4, s) = (islast(orb, s) ? nothing : next(orb, s))
@inline Base.length(orb::Orbit4) = length(orb.elts)

@inline function Base.:(==)(o1::Orbit4, o2::Orbit4)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit4) = orb.pt
