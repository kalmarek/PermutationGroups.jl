using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using LinearAlgebra
using Test
using Revise
using PermutationGroups

include("src/sff.jl")

G = SymmetricGroup(4)
S = gens(G)
ccG = conjugacy_classes(G)

# Multiplication tables for conjugacy classes.
Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]
# @show Ns

F = GF(PermutationGroups.dixon_prime(ccG))

M = matrix(F, Ns[2])
R, T = PolynomialRing(F, string(gensym()))
p = minpoly(R, M)
factors = square_free_factorization(p)

