#########################################################
# Elementary action/properties of permutation actions

@inline @inbounds function Base.:^(n::Integer, p::Perm)
    if 1 <= n <= length(p.d)
        return oftype(n, p.d[n])
    end
    return n
end

@inline Base.:^(v::Tuple, p::Perm) = ntuple(i->v[i^p], length(v))
@inline function Base.:^(v::AbstractVector, p::Perm)
    res = similar(v);
    @inbounds for i in eachindex(v)
        res[i] = v[i^p]
    end
    return res
end

@inline fixes(p::GroupElem, pt, op=^) = op(pt, p) == pt
@inline fixes(p::Perm, v::AbstractVector, op=^) =
    all( i-> v[i] == v[op(i, p)], eachindex(v))

@inline Base.isone(p::Perm) = all(i->i^p == i, 1:degree(p))

fixedpoints(p::Perm, range=1:degree(p), op=^) = [i for i in range if fixes(p, i, op)]
nfixedpoints(p::Perm, range=1:degree(p), op=^) = count(i->fixes(p, i, op), range)

for (fname, findname) in [(:firstmoved, :findfirst), (:lastmoved, :findlast)]
    @eval begin
        function $fname(p::Generic.Perm{I}, op=^) where I
            k = $findname(i -> i != p.d[i], eachindex(p.d))
            k === nothing && return zero(I)
            return I(k)
        end
    end
end

function Generic.emb(p::Generic.Perm{I}, n) where I
    return Generic.emb!(Perm(I(n)), p, 1:degree(p))
end

"""
    conj!(out::Perm, h::Perm, g::Perm)
Computes the conjugation action of `g` on `h` and stores the result in `out`.

The action is understood to be `h → g^-1*h*g`. `out` will be unaliased, if necessary.
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

degree(p::Perm{I}) where I = I(length(p.d))
