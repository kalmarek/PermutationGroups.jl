@testset "PermGroups" begin
    @test PermGroup([perm"(1)"]) isa PermGroup
    @test PermGroup(perm"(1)") isa PermGroup
    G = PermGroup(perm"(1)")
    @test isdefined(G, :stabchain) == false
    @test order(G) == 1
    @test order(G) isa BigInt
    @test order(Int, G) isa Int
    @test isdefined(G, :stabchain) == true

    G = PermGroup(perm"(1,2,3,4)", perm"(1,2)(4)")
    @test StabilizerChain(G) isa StabilizerChain
    sc = StabilizerChain(G)
    @test length(sc) == 3
    @test order(sc) == factorial(4)
    @test order(sc) isa BigInt
    @test order(Int, sc) isa Int

    @test length(PG.basis(G)) == 3
    @test order(G) == factorial(4)

    SN(n) = [Perm(circshift(collect(1:n), -1)), Perm([[2, 1]; 3:n])]

    a = PG.Perms.perm(gens(G, 1))
    @test gens(G, 1) * a == a * gens(G, 1) == a^2

    for N in 6:20
        S = SN(N)
        S = [S; [Perm(randperm(N)) for _ in 1:3]]
        G = PermGroup(S)
        @test order(G) == factorial(N)
    end

    G = PermGroup(perm"(1,2,3,4)", perm"(1,2)")

    m = match(r"order (\d+)", string(G))
    @test m === nothing

    @test perm"(1,3)" in G
    @test perm"(1,5)" ∉ G

    m = match(r"order (\d+)", string(G))
    @test parse(Int, m.captures[1]) == 24

    A = PermGroup(perm"(1,2,3)", perm"(2,3,4)")
    @test order(A) == 12
    m = match(r"order (\d+)", string(A))
    @test parse(Int, m.captures[1]) == 12

    @test perm"(1,2)" ∉ A
    @test perm"(1,3,4)" in A

    @testset "Iterate over PermGroup" begin
        K1 = PermGroup(Perm{UInt32}[perm"(5,6)", perm"(1,2,3,4,5,6)"]) # Symmetric group on 6 symbols
        elements = [g for g in K1]
        @test elements isa Vector{<:Permutation{Perm{UInt32}}}
        uniq_elements = unique(elements)
        @test order(K1) == length(uniq_elements) == 720
        @test uniq_elements == elements

        K2 = PermGroup(Perm{UInt16}[perm"(3,4,5)", perm"(1,2,3,4,5)"]) # Alternating group on 5 symbols
        elements = [g for g in K2]
        @test elements isa Vector{<:Permutation{Perm{UInt16}}}
        uniq_elements = unique(elements)
        @test order(K2) == length(uniq_elements) == 60
        @test uniq_elements == elements
    end
end
