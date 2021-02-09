@testset "StabChain unit tests" begin
    a,b = [perm"(1,2,4,3)(5)", perm"(1,2,5,4)"]
    sc = StabilizerChain([a,b])
    @test length(sc) == 1
    pt, S, Δ = sc[1]
    @test pt == 1
    @test S isa Vector{<:Perm}
    @test S == [a,b]
    @test Δ isa Schreier
    @test collect(Δ) == [1,2,4,5,3]

    # first point on the orbit:
    u = Δ[1]
    # the first Schreier generator is trivial
    @test isone(u*a*getinv(Δ, 1^a))
    # but the second is not
    @test u*b*getinv(Δ, 1^b) == perm"(2,5)(3,4)"
    c = perm"(2,5)(3,4)"
    h, new_depth = sift(c, sc)
    @test h == c
    @test new_depth == 2

    # so we extend the base
    @test firstmoved(c) == 2
    push!(sc, firstmoved(c))

    # pushing extends base...
    @test length(sc) == 2
    @test length(sc.base) == 2
    @test sc.base[2] == firstmoved(c)
    @test length(sc.sgs) == 2

    # and allocates generating set and the transversal
    @test sc.sgs[2] == Perm{Int}[] # we've just extended basis here!
    @test length(sc.transversals) == 2 # the sc.transversals[2] is right now a garbage

    # add c to generators at depth new_depth=2
    push!(sc, c, new_depth, recompute=true)
    # pushing the generator adds to the basis ath the given depth
    @test sc.sgs[2] == [c]
    # and (re)computes the orbit at the given depth:
    @test sc.transversals[2] == Schreier([c], sc.base[2])

    # c should sift now to the identity at depth 3
    x, d = sift(c, sc)
    @test isone(x)
    @test d == 3

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
    @test sift(d, sc) == (d, 2)

    push!(sc, d, 2, recompute=true)

    @test collect(sc.transversals[2]) == [2,5,3,4]

    x, d = sift(a^2*b, sc)
    @test isone(x)
    @test d == 3

    @test order(sc) == 20

    m = match(r"size (\d+)", string(sc))
    @test parse(Int, m.captures[1]) == 20

    m = match(r"Orbit:\s+(\[.*\])", sprint(show, MIME"text/plain"(), sc))
    @test Meta.parse(m.captures[1]) |> eval == collect(sc.transversals[1])
end

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

    @test length(base(G)) == 3
    @test order(G) == factorial(4)

    SN(n) = [Perm(circshift(collect(1:n), -1)), Perm([[2,1];3:n])]

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
end

@testset "Iterate over PermGroup" begin
	K1 = PermGroup([perm"(5,6)", perm"(1,2,3,4,5,6)"]) # Symmetric group on 6 symbols
	elements = [g for g in K1]
    @test elements isa Vector{Perm{Int64}}
	uniq_elements = unique(elements)
    @test order(K1) == length(uniq_elements) == 720
    @test uniq_elements == elements

	K2 = PermGroup(Perm{Int16}[perm"(3,4,5)", perm"(1,2,3,4,5)"]) # Alternating group on 5 symbols
	elements = [g for g in K2]
    @test elements isa Vector{Perm{Int16}}
	uniq_elements = unique(elements)
    @test order(K2) == length(uniq_elements) == 60
	@test uniq_elements == elements
end
