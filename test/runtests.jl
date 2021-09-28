using PermutationGroups
using GroupsCore
using Random
using Test
using BenchmarkTools

import AbstractAlgebra
SymmetricGroup = AbstractAlgebra.SymmetricGroup

@testset "PermutationGroups" begin
    include("orbit.jl")
    include("schreier.jl")
    include("stabchain.jl")
    include("conj.jl")
    include("groups_interface.jl")

    include("benchmark.jl")
end
