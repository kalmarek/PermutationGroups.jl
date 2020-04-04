using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using LinearAlgebra
using Test
using Revise
using PermutationGroups


chars = let G = SymmetricGroup(4)
    S = gens(G)
    ccG = conjugacy_classes(G)

    ccG = [
        Orbit(S, one(G)),
        Orbit(S, perm"(1,2)(3,4)"),
        Orbit(S, perm"(1,2)(4)"),
        Orbit(S, perm"(1,2,3)(4)"),
        Orbit(S, perm"(1,2,3,4)"),
    ]

    @assert sum(length, ccG) == order(G)

    # Multiplication tables for conjugacy classes.
    Ns = (PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG))

    F = GF(PermutationGroups.dixon_prime(ccG))
    # F = GF(41)
    eig_decomposition = PermutationGroups.sd_basis(Ns, F)
    display(eig_decomposition)
    for N in Ns
        # basis actually diagonalizes all Ns
        m = PermutationGroups._change_basis(matrix(F,N), eig_decomposition.basis)
        display(m)
        @test isdiag(Matrix(m))
    end

    chars = [PermutationGroups.Character(
        vec(eig_space), ccG) for eig_space in eig_decomposition]

    # orthogonality of characters:
    @test [dot(χ, ψ) for χ in chars, ψ in chars] ==
        Matrix{Int}(I, length(chars), length(chars))

    chars
end

hcat([χ.vals for χ in chars])


chars = let G = PermGroup([perm"(1,2,4,5,3)", perm"(2,5,3,4)"]);
    @test order(G) == 20

    ccG = conjugacy_classes(G)
    @test sum(length, ccG) == 20

    S = gens(G)
    a,b = S

    ccG = [
        Orbit(S, one(G)),
        Orbit(S, a),
        Orbit(S, b),
        Orbit(S, b^2),
        Orbit(S, b^3),
        ]

    @test sum(length, ccG) == 20
    Ns = (PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG))

    F = GF(PermutationGroups.dixon_prime(ccG))
    # F = GF(41)
    eig_decomposition = PermutationGroups.sd_basis(Ns, F)
    display(eig_decomposition)
    for N in Ns
        # basis actually diagonalizes all Ns
        m = PermutationGroups._change_basis(matrix(F,N), eig_decomposition.basis)
        display(m)
        @test isdiag(Matrix(m))
    end

    chars = [PermutationGroups.Character(
        vec(eig_space), ccG) for eig_space in eig_decomposition]

    # orthogonality of characters:
    @test [dot(χ, ψ) for χ in chars, ψ in chars] ==
        Matrix{Int}(I, length(chars), length(chars))

    chars
end

hcat([χ.vals for χ in chars])

end






# GAP output, read row-wise
# [ [      1,      1,      1,      1,      1 ],
# [      1,      1,     -1,     -1,      1 ],
# [      1,     -1,  -E(4),   E(4),      1 ],
# [      1,     -1,   E(4),  -E(4),      1 ],
# [      4,      0,      0,      0,     -1 ] ]
