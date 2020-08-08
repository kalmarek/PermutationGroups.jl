@testset "Conjugacy classes" begin
    G = PermutationGroups.AbstractAlgebra.SymmetricGroup(4)
    S = [perm"(1,2)(4)", perm"(1,2,3,4)"]
    pt = perm"(4)"
    @test Orbit(S, pt) isa Orbit
    @test length(Orbit(S, pt)) == 1

    pt = perm"(1,2)(4)"
    @test Orbit(S, pt) isa Orbit
    @test length(Orbit(S, pt)) == 6
    @test all(==(permtype(pt)), permtype.(Orbit(S, pt)))

    pt = perm"(1,2)(4)"
    @test Schreier(S, pt) isa Schreier
    tr = Schreier(S, pt)
    @test all(pt^tr[p] == p for p in tr)

    pt = perm"(2,3,4)"
    @test Schreier(S, pt) isa Schreier
    tr = Schreier(S, pt)
    @test all(pt^tr[p] == p for p in tr)

    @test gens(G) == S
    @test conjugacy_classes(G) isa Vector{<:AbstractOrbit}

    ccs = conjugacy_classes(G)
    for cc in ccs
        @test all(permtype(g) == permtype(first(cc)) for g in cc)
    end
end
