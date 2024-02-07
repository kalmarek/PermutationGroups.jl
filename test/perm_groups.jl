@testset "PermGroups" begin
    @test PermGroup([perm"(1)"]) isa PermGroup
    @test PermGroup(perm"(1)") isa PermGroup
    G = PermGroup(perm"(1)")
    @test isdefined(G, :stabchain) == false
    @test order(G) == 1
    @test order(G) isa BigInt
    @test order(Int, G) isa Int
    @test isdefined(G, :stabchain) == true

    @test_throws AssertionError PermGroup()

    G = PermGroup(perm"(1,2)", perm"(1,2,3,4)")
    @test StabilizerChain(G) isa StabilizerChain
    sc = StabilizerChain(G)
    @test length(sc) == 3
    @test order(sc) == factorial(4)
    @test order(sc) isa BigInt
    @test order(Int, sc) isa Int

    @test length(PG.basis(G)) == 3
    @test order(G) == factorial(4)
    @test AP.degree(G) == 4

    H = PermGroup(Permutation(perm"(1,2,3)", G))
    @test order(Int, H) == 3
    @test AP.degree(H) == 3

    SN(n) = [Perm([2:n; 1]), Perm([[2, 1]; 3:n])]

    a = AP.perm(gens(G, 1))
    @test gens(G, 1) * a == a * gens(G, 1) == a^2

    for N in 6:20
        S = SN(N)
        S = [S; [Perm(randperm(N)) for _ in 1:3]]
        G = PermGroup(S)
        @test order(G) == factorial(N)
    end

    G = PermGroup(perm"(1,2)", perm"(1,2,3,4)")
    showstr(x) = sprint((io, x) -> show(IOContext(io, :limit => true), x), x)
    replstr(x) = sprint(
        (io, x) ->
            show(IOContext(io, :limit => true), MIME("text/plain"), x),
        x,
    )
    @test showstr(G) == "PermGroup( (1,2), (1,2,3,4) )"
    @test replstr(G) ==
          "Permutation group on 2 generators generated by\n (1,2)\n (1,2,3,4)"

    m = match(r"order (\d+)", sprint(show, MIME"text/plain"(), G))
    @test m === nothing

    @test perm"(1,3)" in G
    @test perm"(1,5)" ∉ G

    order(Int, G)

    m = match(r"order (\d+)", sprint(show, MIME"text/plain"(), G))
    @test parse(Int, m.captures[1]) == 24

    A = PermGroup(perm"(1,2,3)", perm"(2,3,4)")
    @test order(A) == 12
    m = match(r"order (\d+)", sprint(show, MIME"text/plain"(), A))
    @test parse(Int, m.captures[1]) == 12

    @test perm"(1,2)" ∉ A
    @test perm"(1,3,4)" in A

    @testset "Schreier-Sims: add generator even when orbit is unchanged" begin
        H = PermGroup([perm"(1,5,4,3,2)", perm"(1,4,5,3)"])
        # when missing a generator on the second level order is 90
        @test order(Int, H) == 120
        @test perm"(2,3)" in H
    end

    @testset "power_by_cycles" begin
        G = PermGroup(Perm([2:32; 1]))
        g = gens(G, 1)
        @test isone(g^32)
        @test (g^8)^3 == g^24
        @test PermutationGroups.AP.power_by_cycles(g, 8) ==
              Base.power_by_squaring(g, 8)
    end

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

    G = PermGroup(perm"(1,4,6)(3,5)", perm"(1,5,4,3)")
    @test order(Int, G) == 120

    @testset "Random subgroups of Sym(n): conjugacy classes" begin
        function conjugacy_classes_orbit(grp::GroupsCore.Group)
            id = one(grp)
            S = gens(grp)
            ordG = order(Int, grp)

            cclasses = [PermutationGroups.Orbit([id])]
            elts_counted = 1

            for g in grp
                any(ccl -> g ∈ ccl, cclasses) && continue
                ccl_g = PermutationGroups.Orbit(g, S, conj)
                elts_counted += length(ccl_g)
                push!(cclasses, ccl_g)
                elts_counted == ordG && break
            end

            elts_counted == ordG || @warn "$elts_counted ≠ $ordG"
            return cclasses
        end

        for n in 2:6
            G = PermGroup(perm"(1,2)", Perm([2:n; 1]))
            for _ in 1:10
                S = rand(G, 2)
                H = PermGroup(S)
                ccls = conjugacy_classes_orbit(H)
                @test sum(length, ccls) == order(Int, H) ≤ factorial(n)
            end
        end
    end
end
