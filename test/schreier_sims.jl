@testset "Schreier-Sims step by step" begin
    PG = PermutationGroups
    S = [perm"(1,3,5,7)(2,4,6,8)", perm"(1,3,8)(4,5,7)"]
    P = eltype(S)
    Tr = PG.Transversal(P)

    sc = PG.StabilizerChain{P,Tr}()
    @test PG.istrivial(sc)

    g = sift(sc, S[1])
    @test g == S[1] # since sc is trivial
    # so we need to extend the chain:
    push!(sc.gens, g)
    @test firstmoved(g) == 1
    sc.transversal = Tr(firstmoved(g), g)
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
    sc.stabilizer.transversal = Tr(firstmoved(schr_gen), schr_gen)
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
end
