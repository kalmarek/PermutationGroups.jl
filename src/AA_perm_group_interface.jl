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

Base.deepcopy_internal(g::Perm, stackdict::IdDict) =
    Perm(Base.deepcopy_internal(g.d, stackdict), false)
