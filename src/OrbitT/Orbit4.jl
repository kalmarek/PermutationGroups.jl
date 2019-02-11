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
    orb.elts[orb.last] = (orb.elts[orb.last][1], new)
    orb.last = new
    push!(orb.elts, new => (val, nothing))
    return orb
end

@inline Base.push!(orb::Orbit4{T,Nothing}, n::T) where T = push!(orb, (n,nothing))

# o is a point in orbit, return (value, nextptr)
@inline Base.getindex(orb::Orbit4, o) = orb.elts[o][1]

@inline Base.iterate(orb::Orbit4) = orb.pt, orb.pt
@inline function Base.iterate(orb::Orbit4{T, S}, s) where {T,S} # s is the previous elt
    if s == orb.last
        return nothing
    end
    x = orb.elts[s][2]::T
    return (x, x)
end

@inline Base.iterate(orb::Orbit4, ::Nothing) = nothing
@inline Base.length(orb::Orbit4) = length(orb.elts)
@inline function Base.:(==)(o1::Orbit4, o2::Orbit4)
    o1.pt == o2.pt && o1.last == o2.last || return false
    return o1.elts == o2.elts
end

@inline Base.first(orb::Orbit4) = orb.pt
