@testset "StabChain unit tests" begin
    a, b = [perm"(1,2,4,3)(5)", perm"(1,2,5,4)"]
    P_t = typeof(a)
    Tr_t = Transversal(P_t)

    sc = StabilizerChain{P_t,Tr_t}()

    push!(sc.gens, a, b)
    sc.stabilizer = StabilizerChain{P_t,Tr_t}()
    PG.recompute_transversal!(sc)

    @test length(sc) == 1
    PG.point(sc) == 1
    @test PG.gens(sc) == [a, b]
    Δ = PG.transversal(sc)
    @test Δ isa Transversal
    @test collect(Δ) == [1, 2, 4, 5, 3]

    # first point on the orbit:
    u = Δ[1]
    # the first Schreier generator is trivial
    @test isone(u * a * inv(Δ[1^a]))
    # but the second is not
    @test u * b * inv(Δ[1^b]) == perm"(2,5)(3,4)"
    c = u * b * inv(Δ[1^b])
    h = PG.sift(sc, c)
    @test h == c

    # so we extend the base
    @test AP.firstmoved(c, Base.OneTo(AP.degree(c))) == 2
    # add c to generators
    push!(PG.gens(sc.stabilizer), c)
    # (re)computes the orbit of the stabilizer
    PG.recompute_transversal!(PG.stabilizer(sc))
    sc.stabilizer.stabilizer = StabilizerChain{P_t,Tr_t}()

    # pushing extends base...
    @test length(sc) == 2
    @test length(PG.basis(sc)) == 2
    @test PG.basis(sc)[2] == AP.firstmoved(c, Base.OneTo(AP.degree(c)))

    @test PG.sgs(sc) == [a, b, c]

    # c should sift now to the identity at depth 3
    @test isone(sift(sc, c))
    @test isone(sift(sc, a * b))
    @test isone(sift(sc, c * a))

    u = Δ[2]
    @test u == a
    # Schreier generators are trivial
    @test u * a == Δ[2^a]
    @test u * b == Δ[2^b]
    # so we move to the next point on the orbit
    u = Δ[4] # β = 4
    @test u * a == Δ[4^a]
    @test u * b * inv(Δ[4^b]) == a^2 * b
    d = a^2 * b
    @test sift(sc, d) == d

    # since d stabilizes the whole orbit
    @test Set(collect(PG.orbit(sc)) .^ (d)) == Set(collect(PG.orbit(sc)))
    # we add it directly to the stabilizer of sc
    push!(PG.gens(PG.stabilizer(sc)), d)
    PG.recompute_transversal!(PG.stabilizer(sc))

    @test collect(PG.transversal(PG.stabilizer(sc))) == [2, 5, 3, 4]

    @test isone(sift(sc, a^2 * b))
    @test order(sc) == 20

    sc_str = sprint(show, sc)

    m = match(r"order (\d+)", sc_str)
    @test parse(Int, m.captures[1]) == 20

    m = match(r"Orbit{(.*)}: (\[.*\])", sprint(show, MIME"text/plain"(), sc))
    @test Meta.parse(m.captures[2]) |> eval == collect(PG.transversal(sc))

    @testset "Iteration over leafs of stabilizer chain tree" begin
        a = perm"(1,2,4)(5,3,6)"

        P_t = typeof(a)
        Tr_t = SchreierTransversal(P_t)

        sc = StabilizerChain{P_t,Tr_t}()

        @test length(sc) == 0
        @test length(PG.leafs(sc)) == 1
        @test collect(PG.leafs(sc)) == [one(P_t)]

        push!(sc.gens, a)
        sc.stabilizer = StabilizerChain{P_t,Tr_t}()
        PG.recompute_transversal!(sc)

        @test length(sc) == 1
        @test length(PG.leafs(sc)) == order(Int, a)
        @test collect(PG.leafs(sc)) == [a^i for i in 0:order(Int, a)-1]
    end
end
