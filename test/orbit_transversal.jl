@testset "Orbits, Transversals" begin
    @test PermutationGroups.Orbit([1, 2, 3]) isa Orbit
    gens = [perm"(4,7,5,6)", perm"(1,9,13,8,6,5,7,4)(10,14,12,11)"]
    @test gens isa Vector{<:Perm{UInt16}}

    pt = 4
    orb = PermutationGroups.Orbit(pt, gens)
    @test orb isa Orbit{Int64}
    @test Orbit(UInt32(pt), gens) isa Orbit{UInt32}
    @test Orbit{Int16}(pt, gens) isa Orbit{Int16}
    @test pt in orb
    @test !(10 in orb)
    @test length(orb) == 8
    @test first(orb) == pt

    tr = Transversal(pt, gens)
    @test length(tr) == length(orb)
    @test tr isa Transversal{Int,eltype(gens)}

    @test Transversal(UInt32(pt), gens) isa Transversal{UInt32,eltype(gens)}
    @test Transversal{Int16}(pt, gens) isa Transversal{Int16}

    @test_throws PG.NotInOrbit tr[10]

    schr = SchreierTransversal(pt, gens)
    @test length(schr) == length(orb)
    @test schr isa SchreierTransversal{Int,eltype(gens)}

    @test SchreierTransversal(UInt32(pt), gens) isa
          SchreierTransversal{UInt32,eltype(gens)}
    @test SchreierTransversal{Int16}(pt, gens) isa SchreierTransversal{Int16}

    @test_throws PG.NotInOrbit tr[10]

    @testset "Schreier Vectors and Stabilizers" begin
        Random.seed!(1)
        SIZE = 30
        S = [Perm{UInt16}(randperm(30)) for _ in 1:3]

        pt = 1
        tr = Transversal(pt, S)
        @test tr[pt] == one(S[1])

        @test PermutationGroups.coset_representative(pt, tr) == one(S[1])
        g, h, k = S
        @test PermutationGroups.coset_representative(pt^g, tr) == g
        @test PermutationGroups.coset_representative((pt^g)^h, tr) == g * h
        @test PermutationGroups.coset_representative(pt^(g * h), tr) == g * h

        schr = SchreierTransversal(pt, S)
        @test schr[pt] == one(S[1])

        @test PermutationGroups.coset_representative(pt^g, schr) == g
        @test PermutationGroups.coset_representative((pt^g)^h, schr) == g * h
        @test PermutationGroups.coset_representative(pt^(g * h), schr) == g * h

        # @test Word(inv.(S), schr.orb, pt) == Word(Int[])
        # @test Word(inv.(S), schr.orb, pt^g) == Word(Int[1])
        # @test Word(inv.(S), schr.orb, (pt^g)^h) == Word(Int[2, 1])
        # @test Word(inv.(S), schr.orb, pt^(g * h)) == Word(Int[2, 1])

        z = g * h * k
        δ = pt^z
        # @test Word(inv.(S), schr.orb, δ, ^) == Word([3, 2, 1])

        @test schr[δ] == z
        @test all([pt^schr[o] == o for o in schr])
        @test all(schr[o] == tr[o] for o in tr)
    end

    @testset "allocations" begin
        gens = [perm"(4,7,5,6)", perm"(1,9,13,8,6,5,7,4)(10,14,12,11)"]
        init_pt = 4
        tr = Transversal(init_pt, gens)

        for pt in tr
            g = tr[pt]
            k = @allocated tr[pt]
            @test k == 0
            @test init_pt^g == pt
        end

        schtr = SchreierTransversal(init_pt, gens)
        for pt in schtr
            g = schtr[pt]
            k = @allocated schtr[pt]
            if pt == first(schtr)
                @test k == 0
            elseif v"1.10" ≤ VERSION < v"1.11"
                if pt == last(schtr)
                    @test k == 432
                else
                    @test k ∈ (128, 144)
                end
            end
            @test init_pt^g == pt
        end
    end
end
