using Test
using GroupsCore
using PermutationGroups
using Random

const PG = PermutationGroups
import AbstractPermutations

include(joinpath(pkgdir(AbstractPermutations), "test", "abstract_perm_API.jl"))

@testset "PermutationGroups" begin
    @testset "Perms" begin
        abstract_perm_interface_test(Perm)
        abstract_perm_interface_test(SPerm{8})
    end

    include("orbit_transversal.jl")
    include("conj_action.jl")

    include("stabchain.jl")
    include("schreier_sims.jl")

    include("perm_groups.jl")
    include("groups_interface.jl")
end

include("benchmark.jl")
