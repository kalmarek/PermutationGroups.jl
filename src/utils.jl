#########################################################
# Elementary action/properties of permutation actions
#########################################################

@inline @inbounds function Base.:^(n::Integer, p::perm)
    if n <= length(p.d)
        return oftype(n, p.d[n])
    end
    return n
end

@inline Base.:^(v::Tuple, p::perm) = ntuple(i->v[i^p], length(v))
@inline Base.:^(v::Vector, p::perm) = [v[i^p] for i in eachindex(v)]

@inline fixes(p::GroupElem, pt, op=^) = op(pt, p) == pt
@inline fixes(p::perm, v::AbstractVector) = all(v[i] == v[i^p] for i in eachindex(v))

@inline Base.isone(p::perm) = all(i->first(i)==last(i), enumerate(p.d))

fixedpoints(p::perm, range=1:length(p.d)) = [i for i in range if fixes(p, i)]

for (fname, findname) in [(:firstmoved, :findfirst), (:lastmoved, :findlast)]
    @eval begin
        function $fname(p::perm{I}, op=^) where I
            k = $findname(i -> i != p.d[i], eachindex(p.d))
            k == nothing && return zero(I)
            # isnothing(k) && return zero(I)
            return I(k)
        end
    end
end

function gensstring(gens::Vector{<:perm}; width=96)
    str = ""
    ellipsis = " … "

    str = join(gens, ", ")
    if length(str) > width
        str = str[1:width - length(ellipsis)] * ellipsis
    end
    return str
end

AbstractAlgebra.degree(p::perm) = length(p.d)
function Generic.emb(p::perm{I}, n) where I
    return Generic.emb!(perm(I(n)), p, 1:degree(p))
end