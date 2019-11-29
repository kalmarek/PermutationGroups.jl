#########################################################
# Elementary action/properties of permutation actions
#########################################################

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

AbstractAlgebra.degree(p::Perm{I}) where I = I(length(p.d))
function Generic.emb(p::Generic.Perm{I}, n) where I
    return Generic.emb!(Perm(I(n)), p, 1:degree(p))
end

# TODO: move to AbstractAlgebra
function fastmul!(out::Perm, g::Perm, h::Perm)
   out = (out === h ? similar(out) : out)
   @inbounds for i in eachindex(out.d)
      out[i] = h[g[i]]
   end
   return out
end

Base.one(G::Generic.PermGroup) = G()
