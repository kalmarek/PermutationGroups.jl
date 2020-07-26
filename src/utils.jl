#########################################################
# Elementary action/properties of permutation actions

@inline @inbounds function Base.:^(n::Integer, p::Perm)
    if n <= length(p.d)
        return oftype(n, p.d[n])
    end
    return n
end

@inline Base.:^(v::Tuple, p::Perm) = ntuple(i->v[i^p], length(v))
@inline Base.:^(v::Vector, p::Perm) = [v[i^p] for i in eachindex(v)]

@inline fixes(p::GroupElem, pt, op=^) = op(pt, p) == pt
@inline fixes(p::Perm, v::AbstractVector) = all(v[i] == v[i^p] for i in eachindex(v))

@inline Base.isone(p::Perm) = all(i->first(i)==last(i), enumerate(p.d))

fixedpoints(p::Perm, range=1:length(p.d)) = [i for i in range if fixes(p, i)]

for (fname, findname) in [(:firstmoved, :findfirst), (:lastmoved, :findlast)]
    @eval begin
        function $fname(p::Generic.Perm{I}, op=^) where I
            k = $findname(i -> i != p.d[i], eachindex(p.d))
            k == nothing && return zero(I)
            # isnothing(k) && return zero(I)
            return I(k)
        end
    end
end

function gensstring(gens::AbstractVector{<:Perm}; width=96)
    str = ""
    ellipsis = " … "

    str = join(gens, ", ")
    if length(str) > width
        str = str[1:width - length(ellipsis)] * ellipsis
    end
    return str
end

#########################################################
# Misc functions that should go to AbstractAlgebra

import Base: one, conj
import AbstractAlgebra: degree

Base.one(G::Generic.SymmetricGroup{I}) where I = Perm(G.n)
Base.one(g::Perm{I}) where I = Perm(degree(g))

AbstractAlgebra.degree(p::Perm{I}) where I = I(length(p.d))

function Generic.emb(p::Generic.Perm{I}, n) where I
    return Generic.emb!(Perm(I(n)), p, 1:degree(p))
end

"""
    conj!(out::Perm, h::Perm, g::Perm)
Computes the conjugation action of `g` on `h` and stores the result in `out`.
The action is understood to be `h → g^-1*h*g`.
`out` will be unaliased, if necessary.
"""
function Base.conj!(out::Perm, h::Perm, g::Perm)
    if out === h
        out = deepcopy(out)
    end
    @inbounds for i in 1:degree(g)
        out[g[i]] = g[h[i]]
    end
    return out
end

Base.conj(h::GroupElem, g::GroupElem) = conj!(h, h, g)
Base.:(^)(h::Perm, g::Perm) = conj(h,g)

function AbstractAlgebra.gens(G::Generic.SymmetricGroup)
    a = one(G)
    b = one(G)
    circshift!(b.d, a.d, -1)
    a.d[1], a.d[2] = 2, 1
    return [a, b]
end

