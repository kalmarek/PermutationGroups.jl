@testset "gens for SymmetricGroup" begin
    @test gens(PermutationGroups.SymmetricGroup(1)) == [Perm(1)]
    @test gens(PermutationGroups.SymmetricGroup(2)) == [perm"(1,2)"]
    @test gens(PermutationGroups.SymmetricGroup(3)) == [perm"(1,2)", perm"(1,2,3)"]
    @test gens(PermutationGroups.SymmetricGroup(4)) == [perm"(1,2)", perm"(1,2,3,4)"]
end
