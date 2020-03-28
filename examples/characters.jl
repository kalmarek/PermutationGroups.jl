using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using Revise
using PermutationGroups

using Test

sd = let G = SymmetricGroup(4)
    S = gens(G)
    ccG = conjugacy_classes(G)

    # Multiplication tables for conjugacy classes.
    Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]

    F = GF(PermutationGroups.dixon_prime(ccG))
    basis = PermutationGroups.sd_basis(Ns, F)
end


sd = let G = PermGroup([perm"(1,2,4,5,3)", perm"(2,5,3,4)"]);
    @test order(G) == 20

    S = gens(G)
    a,b = S

    ccG = [Orbit(S, one(G)),
        Orbit(S, a),
        Orbit(S, b),
        Orbit(S, b^2),
        Orbit(S, b^3),
        ]

    @test sum(length, ccG) == 20
    Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]

    # F = GF(PermutationGroups.dixon_prime(ccG))
    F = GF(101)
    basis = PermutationGroups.sd_basis(Ns, F)
end


A = matrix(GF(691), [
     0   1  0  0
     0   0  5  0
     0   0  0  5
     5   0  0  0]
    )

# GAP output, read row-wise
# [ [      1,      1,      1,      1,      1 ],
  # [      1,      1,     -1,     -1,      1 ],
  # [      1,     -1,  -E(4),   E(4),      1 ],
  # [      1,     -1,   E(4),  -E(4),      1 ],
  # [      4,      0,      0,      0,     -1 ] ]




