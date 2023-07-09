using Test
using GroupsCore
using PermutationGroups
using Random

const PG = PermutationGroups

@testset "Perms" begin
    include("AbstractPerm_interface.jl")
    include("perm_interface.jl")

    abstract_perm_interface_test(APerms.APerm)
    abstract_perm_interface_test(Perm)

    @test_throws AssertionError Perm([1, 2, 3, 1])

    # more implementations come here …

    include("perm_macro.jl")
end

@testset "PermutationGroups" begin
    include("orbit_transversal.jl")
    include("conj_action.jl")

    include("stabchain.jl")
    include("schreier_sims.jl")

    include("perm_groups.jl")
    include("groups_interface.jl")
end

include("benchmark.jl")
