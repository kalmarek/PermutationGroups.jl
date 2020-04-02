using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using LinearAlgebra
using Test
using Revise
using PermutationGroups

sd = let G = SymmetricGroup(4)
    S = gens(G)
    ccG = conjugacy_classes(G)

    # Multiplication tables for conjugacy classes.
    Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]
    # @show Ns

    F = GF(PermutationGroups.dixon_prime(ccG))

    basis = PermutationGroups.sd_basis(Ns, F)
    ib, b = inv(basis.basis), basis.basis
    for N in Ns
        @test isdiag(Matrix(ib*matrix(F,N)*b))
    end
    basis
end

sd.basis*sd.basis'


sd = let
    G = PermGroup([perm"(1,2,4,5,3)", perm"(2,5,3,4)"]);
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
    Ns = (PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG))

    F = GF(PermutationGroups.dixon_prime(ccG))
    # F = GF(41)
    basis = PermutationGroups.sd_basis(Ns, F)

    ib, b = inv(basis.basis), basis.basis
    for N in Ns
        @test isdiag(Matrix(ib*matrix(F,N)*b))
    end
    basis
end






# GAP output, read row-wise
# [ [      1,      1,      1,      1,      1 ],
# [      1,      1,     -1,     -1,      1 ],
# [      1,     -1,  -E(4),   E(4),      1 ],
# [      1,     -1,   E(4),  -E(4),      1 ],
# [      4,      0,      0,      0,     -1 ] ]
