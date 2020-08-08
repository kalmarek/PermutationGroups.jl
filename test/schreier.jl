@testset "Schreier Vectors and Stabilizers" begin
    Random.seed!(1);
    SIZE=30
    G = PermutationGroups.AbstractAlgebra.SymmetricGroup(SIZE);
    S = [rand(G) for _ in 1:3];

    pt = 1
    tr = Orbit1(Transversal, S, pt)

    import PermutationGroups.Word

    for OrbT in OrbitTypes
        schr = Schreier(OrbT, S, pt, ^)
        @test schr[pt] == one(G)

        g,h,k = S

        @test representative(S, schr.orb, pt) == one(G)
        @test representative(S, schr.orb, pt^g) == g
        @test representative(S, schr.orb, (pt^g)^h) == g*h
        @test representative(S, schr.orb, pt^(g*h)) == g*h

        @test Word(inv.(S), schr.orb, pt) == Word(Int[])
        @test Word(inv.(S), schr.orb, pt^g) == Word(Int[1])
        @test Word(inv.(S), schr.orb, (pt^g)^h) == Word(Int[2,1])
        @test Word(inv.(S), schr.orb, pt^(g*h)) == Word(Int[2,1])

        z = g*h*k
        δ = pt^z
        @test Word(inv.(S), schr.orb, δ, ^) == Word([3,2,1])

        @test schr[δ] == z
        @test all([pt^schr[o] == o for o in schr.orb])
        @test all(schr[o] == tr[o] for o in tr)
    end
end


@testset "Schreier-Sims unit tests" begin
    S = [perm"(1,2,4,3)(5)", perm"(2,5,4)"]
    @test PermutationGroups.initial_bsgs(S)[1] == [1,2]
    @test PermutationGroups.initial_bsgs(S)[2] == [S, [S[2]]]

    a,b = [perm"(1,2,4,3)(5)", perm"(1,2,5,4)"]
    B, S = PermutationGroups.initial_bsgs([a,b])
    @test B == [1]
    @test S[1] == [a,b]

    trs = [Schreier(gens, pt) for (gens, pt) in zip(S, B)]
    @test collect(trs[1]) == [1,2,4,5,3]
    depth = 1
    Δ = trs[depth]

    # first point on the orbit is β = 1:
    u = Δ[1]
    # the first Schreier generator is trivial
    @test isone(u*a*getinv(Δ, 1^a))
    # but the second is not
    @test u*b*getinv(Δ, 1^b) == perm"(2,5)(3,4)"
    c = perm"(2,5)(3,4)"
    @test sift(c, B, trs) == (c, 2)

    # so we extend the basis
    push!(B, firstmoved(c))
    # add c to generators at depth 2
    push!(S, [c])
    # and compute the orbit of the added point at depth 2
    push!(trs, Schreier(S[depth+1], 2))
    depth = 2
    @test sift(c, B, trs) == (perm"(5)", 3)

    u = Δ[2]
    @test u == a
    # Schreier generators are trivial
    @test u*a == Δ[2^a]
    @test u*b == Δ[2^b]
    # so we move to the next point on the orbit
    u = Δ[4] # β = 4
    @test u*a == Δ[4^a]
    @test u*b*getinv(Δ, 4^b) == a^2*b
    d = a^2*b
    @test sift(d, B, trs) == (d, 2)

    push!(S[2], d)
    trs[2] = Schreier(S[2], B[2])
    @test collect(trs[2]) == [2,5,3,4]

    @test sift(a^2*b, B, trs) == (perm"(5)", 3)

    @test prod(length, trs) == 20

    B_new, S_new, trs_new = schreier_sims([a,b])
    @test B_new == B
    @test S_new == S
    @test trs_new == trs

end
