using AbstractAlgebra
using GroupRings

struct CCMatrix{T, C} <: AbstractMatrix{T} # M_r
    cc::Vector{C} # vector of conjugacy classes to fix the order
    r::Int # the index of conjugacy class
    m::Matrix{T} # cache of class coefficients

    function CCMatrix(cc::A, r::Int, T::Type=Int) where {C, 
A<:AbstractVector{C}}
        M = -ones(T, length(cc), length(cc))
        new{T, C}(cc, r, M)
    end
end

Base.size(M::CCMatrix) = size(M.m)
Base.IndexStyle(::Type{<:CCMatrix}) = IndexCartesian()

function Base.getindex(M::CCMatrix{T, C}, s::Integer, t::Integer) where {T, 
        C <: GroupRingElem}
    if isone(-M.m[s,t])
        r = M.r
        g = first(supp(M.cc[t])) # it doesn't matter which we take
        M.m[s,t] = (M.cc[r]*M.cc[s])[g] # c_RST = r*s = t
    # we compute much more above: no we obtain the whole row (? or column?)
    end
    return M.m[s,t]
end

N = 4
G = PermGroup(N)
RG = GroupRing(G)

ccG = let ptypes = [p.part for p in AllParts(N)]
    perms_typed = Dict(g => permtype(g) for g in G)
    ccG = [Set(g for (g,t) in perms_typed if t==pt) for pt in ptypes]
    [sum(RG.(cc)) for cc in ccG]
    # conjugacy classes as all-one elements of GroupRing
end

M = [CCMatrix(ccG, i) for i in 1:length(ccG)];
@assert M[1] == one(M[1])
          
# <a, b> = sum(a(g^-1)*b(g) for g in RepresentativesConjugacyClasses(G) )
