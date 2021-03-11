using PermutationGroups
using Random
using Test
using BenchmarkTools

@testset "PermutationGroups" begin
    include("utils.jl")
    include("orbit.jl")
    include("schreier.jl")
    include("stabchain.jl")
    include("conj.jl")
    
    include("benchmark.jl")
end
