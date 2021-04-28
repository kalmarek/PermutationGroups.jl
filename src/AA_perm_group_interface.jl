# to be deleted when/if AA satisfies Groups Interface

GroupsCore.order(
    ::Type{I},
    G::SymmetricGroup,
) where {I<:Integer} = I(factorial(G.n))

# disambiguation
GroupsCore.order(
    ::Type{I},
    g::Perm,
) where {I<:Integer} =
    I(foldl(lcm, length(c) for c in cycles(g)))

# new methods

function GroupsCore.gens(G::SymmetricGroup{I}) where {I}
    degree(G) == 1 && return [one(G)]
    a, b = one(G), one(G)
    circshift!(a.d, b.d, -1)
    degree(G) == 2 && return [a]
    b.d[1], b.d[2] = 2, 1
    return [a, b]
end

Base.deepcopy_internal(g::Perm, stackdict::IdDict) =
    Perm(Base.deepcopy_internal(g.d, stackdict), false)
