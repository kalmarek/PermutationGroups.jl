@testset "Schreier-Sims step by step" begin
    PG = PermutationGroups
    S = [perm"(1,3,5,7)(2,4,6,8)", perm"(1,3,8)(4,5,7)"]
    P = eltype(S)
    Tr = PG.__schreier_sims_transversal(Transversal, P)

    sc = PG.StabilizerChain{P,Tr}()
    @test PG.istrivial(sc)

    g = sift(sc, S[1])
    @test g == S[1] # since sc is trivial
    # so we need to extend the chain:
    push!(sc.gens, g)
    @test AP.firstmoved(g, Base.OneTo(AP.degree(g))) == 1
    sc.transversal = Tr(AP.firstmoved(g, Base.OneTo(AP.degree(g))), g)
    sc.stabilizer = PG.StabilizerChain{P,Tr}()

    @test !(PG.istrivial(sc))
    @test collect(PG.orbit(sc)) == [1, 3, 5, 7]

    g = sift(sc, S[2])
    @test g == S[2] * inv(sc.transversal[first(sc.transversal)^S[2]])
    # so we need to recalculate the orbit of 1 under [S[1], g]
    push!(sc.gens, g)
    PG.recompute_transversal!(sc)
    @test collect(PG.orbit(sc)) == [1, 3, 5, 6, 7, 8, 4, 2]
    Δ = PG.transversal(sc)
    pt = first(PG.orbit(sc))
    u = Δ[pt]
    @test isone(u)
    @test u * g * inv(Δ[pt^g]) == g # since g fixes pt
    # process the schreier generator
    schr_gen = g
    push!(sc.stabilizer.gens, schr_gen)
    sc.stabilizer.transversal =
        Tr(AP.firstmoved(schr_gen, Base.OneTo(AP.degree(schr_gen))), schr_gen)
    sc.stabilizer.stabilizer = PG.StabilizerChain{P,Tr}()

    @test collect(PG.orbit(sc.stabilizer)) == [2, 8, 7]
    @test PG.order(Int, sc) == 24

    # since we know that order(PermGroup(S)) == 24 we could stop at this moment.

    pt = let orb = PG.orbit(sc)
        _, st = iterate(orb)
        pt, _ = iterate(orb, st)
        pt
    end

    # let's check that the next Schreier generators sift to trivial ones:
    u = Δ[pt]
    for g in PG.gens(sc)
        schr = u * g * inv(Δ[pt^g])
        @test isone(PG.sift(sc, schr))
    end

    for pt in PG.orbit(sc)
        for g in PG.gens(sc)
            schr = u * g * inv(Δ[pt^g])
            @test isone(PG.sift(sc, schr))
        end
    end

    @test_logs (:warn, "schreier_sims using order is not implemented yet") PG.schreier_sims(
        S,
        24,
    )
    sc2 = PG.schreier_sims(S, 24)
    @test PG.order(Int, sc2) == PG.order(Int, sc)
    @test PG.basis(sc2) == PG.basis(sc) == [1, 2]
    @test PG.gens(sc2) == PG.gens(sc)

    @testset "perm from base images" begin
        cube222 = Perm.([
            [1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16],
            [1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16],
            [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 11, 12, 15, 16],
            [16, 8, 14, 6, 12, 4, 10, 2, 15, 7, 13, 5, 11, 3, 9, 1],
            [3, 11, 1, 9, 7, 15, 5, 13, 4, 12, 2, 10, 8, 16, 6, 14],
        ])
        sc = PG.schreier_sims(cube222)
        β = PG.basis(sc)
        @test all(PG.leafs(sc)) do g
            β_im = β .^ g
            return g == AP.perm(sc, β_im)
        end
    end
end
