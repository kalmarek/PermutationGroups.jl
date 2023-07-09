function abstract_perm_interface_test(P::Type{<:PG.AbstractPermutation})
    @testset "AbstractPermutation API test: $P" begin
        @test one(P) isa PermutationGroups.AbstractPermutation
        @test isone(one(P))
        @test 5^one(P) == 5
        @test UInt16(5)^one(P) isa UInt16

        r = [3, 1, 2, 4]
        p = P(r)
        q = P(r[1:3])
        @test degree(p) == 3
        @test degree(q) == 3

        @test p == q
        @test q == Perm(r)

        @test hash(p) == hash(q)
        @test length(unique([p, q])) == 1

        @test collect(PermutationGroups.cycles(p)) == [[1, 3, 2]]
        @test collect(PermutationGroups.cycles(q)) == [[1, 3, 2]]

        @test inv(p) isa PermutationGroups.AbstractPermutation
        @test isone(inv(p) * p)
        @test isone(p * inv(p))

        @test isone(inv(q) * p)
        @test isone(inv(p) * Perm(r))

        a = P([2, 1, 3]) # (1,2)
        b = P([2, 3, 1]) # (1,2,3)

        @test a * b == P([3, 2, 1]) # (1,2)*(1,2,3) == (1,3)
        @test b * a == P([1, 3, 2]) # (1,2,3)*(1,2) == (2,3)
        @test isone(a * a)
        @test isone(b * b * b)

        @test p^q == p
        @test p^Perm(r) == p
        @test p^Perm(r) isa P
        @test Perm(r)^p == p
        @test Perm(r)^p isa Perm

        @test p * q == p^2
        @test p * Perm(r) == p^2
        @test p * Perm(r) * p isa P
        @test p * Perm(r) * p * q isa P

        @test Perm(r) * p isa Perm

        @test (1:5) .^ p == [3, 1, 2, 4, 5]
        @test sprint(show, p) == "(1,3,2)"
    end
end
