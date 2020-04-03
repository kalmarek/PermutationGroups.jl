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
#=
G = SymmetricGroup(4)
S = gens(G)
ccG = conjugacy_classes(G)

# Multiplication tables for conjugacy classes.
Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]
# @show Ns

F = GF(PermutationGroups.dixon_prime(ccG))

for i = 1:length(Ns)
M = matrix(F, Ns[i])
R, T = PolynomialRing(F, string(gensym()))
R, T = PolynomialRing(F, "T")

p = minpoly(R, M)
factors = square_free_factorization(p, T, order(F))
println(factors)
end
=#
