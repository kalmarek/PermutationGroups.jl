using PermutationGroups
using Random
using Test
using BenchmarkTools

@testset "PermutationGroups" begin
    include("orbit.jl")
    include("schreier.jl")
    include("stabchain.jl")
    include("conjugacy_classes.jl")

    include("benchmark.jl")
end
