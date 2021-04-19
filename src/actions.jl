#########################################################
# Elementary actions of Abstract permutations

perm(p::AbstractPerm, op=^) = Perm([op(i, p) for i in Base.OneTo(degree(p))])

@inline fixes(p::GroupElement, pt, op=^) = op(pt, p) == pt
@inline fixes(p::AbstractPerm, v::AbstractVector, op=^) =
    all( i-> v[i] == v[op(i, p)], eachindex(v))

@inline Base.isone(p::AbstractPerm) = all(i->i^p == i, Base.OneTo(degree(p)))

fixedpoints(p::AbstractPerm, range=Base.OneTo(degree(p)), op=^) =
    [i for i in range if fixes(p, i, op)]
nfixedpoints(p::AbstractPerm, range=Base.OneTo(degree(p)), op=^) = count(i->fixes(p, i, op), range)

for (fname, findname) in [(:firstmoved, :findfirst), (:lastmoved, :findlast)]
    @eval begin
        function $fname(p::AbstractPerm, op=^) where I
            k = $findname(i -> i != op(i, p), Base.OneTo(degree(p)))
            k === nothing && return k
            return k
        end
    end
end

@inline Base.:^(v::Tuple, p::AbstractPerm) = ntuple(i->v[i^p], length(v))
@inline function Base.:^(v::AbstractVector, p::AbstractPerm)
    res = similar(v);
    @inbounds for i in eachindex(v)
        res[i] = v[i^p]
    end
    return res
end

# specific definitions

@inline @inbounds function Base.:^(n::Integer, p::Perm)
    if 1 <= n <= length(p.d)
        return oftype(n, p.d[n])
    end
    return n
end

degree(p::Perm{I}) where I = I(length(p.d))
