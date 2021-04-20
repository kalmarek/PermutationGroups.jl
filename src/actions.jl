#########################################################
# Elementary actions of Abstract permutations

Base.eachindex(p::AbstractPerm) = Base.OneTo(degree(p))
perm(p::AbstractPerm, op=^) = Perm([op(i, p) for i in eachindex(p)])

@inline fixes(p::GroupElement, pt, op=^) = op(pt, p) == pt
@inline fixes(p::AbstractPerm, v::AbstractVector, op=^) =
    all( i-> v[i] == v[op(i, p)], eachindex(v))

fixedpoints(p::AbstractPerm, range=eachindex(p), op=^) =
    [i for i in range if fixes(p, i, op)]
nfixedpoints(p::AbstractPerm, range=eachindex(p), op=^) = count(i->fixes(p, i, op), range)

for (fname, findname) in [(:firstmoved, :findfirst), (:lastmoved, :findlast)]
    @eval begin
        function $fname(p::AbstractPerm, op=^)
            k = $findname(i -> i != op(i, p), eachindex(p))
            k === nothing && return k
            return eltype(p)(k)
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

Base.:^(n::Integer, p::Permutation) = n^perm(p)
