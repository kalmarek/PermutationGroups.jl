@testset "misc Permutation functionality" begin
    G = PermGroup([perm"(1,2)", perm"(1,2,3,4)"])
    S = gens(G)

    @test permtype(S[1]) == [2]
    @test permtype(S[2]) == [4]
    @test sign(S[1]) == -1
    @test sign(S[2]) == -1
    @test sign(S[1] * S[2]) == 1
end
