include(joinpath(pathof(GroupsCore), "..", "..", "test", "conformance_test.jl"))

@testset "GroupsCore interface conformance" begin
    G = PermGroup(perm"(1,2,3,4,5)", perm"(1,2)")
    test_Group_interface(G)
    test_GroupElement_interface(rand(G, 2)...)
end
