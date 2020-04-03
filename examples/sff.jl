using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using LinearAlgebra
using Test
using Revise
using PermutationGroups

include("../src/sff.jl")

F = GF(3)
R, T = PolynomialRing(F, "x")
p = T^11 + 2*T^9 + 2*T^8 + T^6 + T^5 + 2*T^3 +2*T^2 + 1
fact = square_free_factorization(p, T, 3)
println(fact)
@test all(x -> p(x) == 0, roots(fact))

G = SymmetricGroup(4)
S = gens(G)
ccG = conjugacy_classes(G)

# Multiplication tables for conjugacy classes.
Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]
F = GF(PermutationGroups.dixon_prime(ccG))

for N in Iterators.rest(Ns, 1)
    M = matrix(F, N)
    @info M
    R, T = PolynomialRing(F, "T")
    p = minpoly(R, M)
    factors = square_free_factorization(p, T, order(F))
    println(roots(factors))
end


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

for N in Iterators.rest(Ns, 1)
    M = matrix(F, N)
    @info M
    R, T = PolynomialRing(F, "T")
    p = minpoly(R, M)
    factors = square_free_factorization(p, T, order(F))
    println(roots(factors))
end


