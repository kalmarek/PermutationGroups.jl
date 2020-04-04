using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using LinearAlgebra
using Test
using Revise
using PermutationGroups


chars = let G = SymmetricGroup(4)
    S = gens(G)
    ccG = conjugacy_classes(G)

    # Multiplication tables for conjugacy classes.
    Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]
    # @show Ns

    F = GF(PermutationGroups.dixon_prime(ccG))

    basis = PermutationGroups.sd_basis(Ns, F)
    ib, b = inv(basis.basis), basis.basis
    for N in Ns
        # basis actually diagonalizes all Ns
        @test isdiag(Matrix(ib*matrix(F,N)*b))
    end
    display(basis)
    chars = [PermutationGroups.Character(vec(v), ccG) for v in basis]

    # orthogonality of characters:
    @test isdiag([dot(χ, ψ) for χ in chars, ψ in chars])

    chars
end

hcat([χ.vals for χ in chars])


chars = let G = PermGroup([perm"(1,2,4,5,3)", perm"(2,5,3,4)"]);
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
    display(basis)
    ib, b = inv(basis.basis), basis.basis
    for N in Ns
        # basis actually diagonalizes all Ns
        @test isdiag(Matrix(ib*matrix(F,N)*b))
    end

    χ = [PermutationGroups.Character(vec(v), ccG) for v in basis]

    # orthogonality of characters:
    @test isdiag([dot(χ[i], χ[j]) for i in 1:length(χ), j in 1:length(χ)])

    χ
end

hcat([χ.vals for χ in chars])

end






# GAP output, read row-wise
# [ [      1,      1,      1,      1,      1 ],
# [      1,      1,     -1,     -1,      1 ],
# [      1,     -1,  -E(4),   E(4),      1 ],
# [      1,     -1,   E(4),  -E(4),      1 ],
# [      4,      0,      0,      0,     -1 ] ]
