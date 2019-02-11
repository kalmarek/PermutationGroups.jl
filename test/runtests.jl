using PermutationGroups
using Random
using AbstractAlgebra
using Test

function manual_orbit(gens::Vector{GEl}, pt::T, op=^) where {GEl <: GroupElem, T}
    G = parent(gens[1])
    Δ = [pt]
    Δ_set = Set(Δ)
    for δ in Δ
        for g in gens
            γ = op(δ, g)
            if !(γ ∈ Δ_set)
                push!(Δ, γ)
                push!(Δ_set, γ)
            end
        end
    end
    return (Δ, Δ_set)
end

function test_orbit(OrbT::Type{T}) where T<: AbstractOrbit
    or = OrbT(1, Symbol("a"))
    @test length(or) == 1
    @test or[1] == :a
    push!(or, (3, :z))
    @test length(or) == 2
    @test collect(or) == [1, 3]
    @test or[1] == :a
    @test or[3] == :z
    push!(or, (2,:c))
    @test length(or) == 3
    push!(or, (5,:b))
    @test length(or) == 4
    @test collect(or) == [1,3,2,5]
    @test or[5] == :b

    orb = OrbT(7)
    @test orb isa OrbT{Int, Nothing}
    @test collect(orb) == [7]
    push!(orb, 5)
    push!(orb, 3)
    @test length(orb) == 3
    @test orb[3] == nothing
    @test orb[5] == nothing
    @test collect(orb) == [7,5,3]
end

@testset "Unit: perms actions" begin
    g = perm"(2,3,5)(8,9)"
    @test fixes(g, 1)
    @test !(fixes(g,2))
    @test fixes(g, 10)
    @test map(i->fixes(g, i), 1:11) == [true, false, false, true, false, true, true, false, false, true, true]
    @test fixedpoints(g) == [1,4,6,7]

    @test firstmoved(g) == 2
    @test firstmoved(perm(5)) == 0
    @test fixedpoints(perm"(1,2,3,4,5)") isa Vector{Int}
    @test fixedpoints(perm(Int8[1,2,3,4,5])) isa Vector{Int}
    @test fixedpoints(perm"(1,3,5)") == [2,4]

    g = perm(Int8[1,3,4,2]) # (1)(2,3,4)
    @test fixes(g, 1) == true
    @test fixes(g, 5) == true
    @test [2,1,1,1]^g isa Vector{Int}

    @test fixes(g, [2,1,1,1]) == true
    @test fixes(g, [2,2,1,1]) == false

    @test firstmoved(g) == 2
    @test firstmoved(g) isa Int8

    # g = (1)(2,3,4)
    @test [1,2,3].^Ref(g) isa Vector{Int}
    @test [1,2,3].^Ref(g) == [1,3,4]
    @test Int16[1,2,3].^Ref(g) isa Vector{Int16}

    @test (1,2,3).^Ref(g) isa NTuple{3, Int}
    @test tuple(Int16[1,2,3]...).^Ref(g) isa NTuple{3,Int16}
    @test (1,2,3).^Ref(g) == (1,3,4)

    @test isone(g) == false
    @test isone(perm([1,2,3]))
end

import PermutationGroups: Orbit1, Orbit2, Orbit3, Orbit4, Orbit5
const OrbitTypes = [Orbit1, Orbit2, Orbit3, Orbit4, Orbit5]

@testset "Orbit Types" begin

    @test Orbit1([1,3,2]) isa Orbit1{Int}

    o = Orbit1([1,3,2])
    @test 1 in o
    @test 4 ∉ o
    push!(o, 4)
    @test 4 in o
    @test collect(o) == [1,3,2,4]

    orb = Orbit2(1, (2,(:b,nothing)), Dict{Int, Tuple{Symbol, Int}}(1=>(:a, 3), 3=>(:c, 2)))
    @test collect(orb) isa Vector{Int}
    @test collect(orb) == [1,3,2]
    @test 1 in orb
    @test 4 ∉ orb
    push!(orb, (4, :z))
    @test 4 in orb
    @test collect(o) == [1,3,2,4]

    orb = Orbit3(1, (2,:b ), Dict{Int, Tuple{Symbol, Int}}(1=>(:a, 3), 3=>(:c, 2)))
    @test collect(orb) isa Vector{Int}
    @test collect(orb) == [1,3,2]
    @test 1 in orb
    @test 4 ∉ orb
    push!(orb, (4, :z))
    @test 4 in orb
    @test collect(o) == [1,3,2,4]

    orb = Orbit4(1, 2, Dict{Int, Tuple{Symbol, Union{Int, Nothing}}}(1=>(:a, 3), 3=>(:c, 2), 2 =>(:b, nothing)))
    @test collect(orb) isa Vector{Int}
    @test collect(orb) == [1,3,2]
    @test 1 in orb
    @test 4 ∉ orb
    push!(orb, (4, :z))
    @test 4 in orb
    @test collect(o) == [1,3,2,4]

    orb = Orbit5(1, 2, Dict{Int, Tuple{Symbol, Int}}(1=>(:a, 3), 3=>(:c, 2), 2 =>(:b, 2)))
    @test collect(orb) isa Vector{Int}
    @test collect(orb) == [1,3,2]
    @test 1 in orb
    @test 4 ∉ orb
    push!(orb, (4, :z))
    @test 4 in orb
    @test collect(o) == [1,3,2,4]
end

@testset "Orbit definitions: elementary ops" begin

    Random.seed!(1);
    SIZE=30
    G = PermGroup(SIZE);
    gens = [rand(G) for _ in 1:3];

    pt = 7
    Δ, _ = manual_orbit(gens, pt, ^)

    for OrbT in OrbitTypes
        test_orbit(OrbT)
        @test collect(OrbT(gens, pt, ^)) == Δ
        t = OrbT(Transversal, gens, pt, ^)
        @test length(t) == SIZE
    end

    t1 = Orbit1(Transversal, gens, pt, ^);

    for OrbT in OrbitTypes
        tr = OrbT(Transversal, gens, pt, ^)
        @test all([tr[o] == t1[o] for o in tr])
        @test all([pt^tr[o] == o for o in tr])
    end
end

@testset "Schreier Vectors and Stabilizers" begin
    Random.seed!(1);
    SIZE=30
    G = PermGroup(SIZE);
    S = [rand(G) for _ in 1:3];

    pt = 1
    tr = Orbit1(Transversal, S, pt)

    import PermutationGroups.Word

    for OrbT in OrbitTypes
        schr = Schreier(OrbT, S, pt, ^)
        @test schr[pt] == G()

        g,h,k = S

        @test representative(S, schr.orb, pt) == G()
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

@testset "StabChain unit tests" begin
    a,b = [perm"(1,2,4,3)(5)", perm"(1,2,5,4)"]
    sc = StabilizerChain([a,b])
    @test length(sc) == 1
    pt, S, Δ = sc[1]
    @test pt == 1
    @test S isa Vector{<:perm}
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
    @test sc.sgs[2] == perm{Int}[] # we've just extended basis here!
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
end

@testset "PrmGroups" begin
    @test PrmGroup([perm"(1)"]) isa PrmGroup
    G = PrmGroup([perm"(1)"])
    @test isdefined(G, :stabchain) == false
    @test order(G) == 1
    @test isdefined(G, :stabchain) == true

    G = PrmGroup([perm"(1,2,3,4)", perm"(1,2)(4)"])
    @test StabilizerChain(G) isa StabilizerChain
    sc = StabilizerChain(G)
    @test length(sc) == 3
    @test order(sc) == factorial(4)

    @test length(base(G)) == 3
    @test order(G) == factorial(4)

    SN(n) = [perm(circshift(collect(1:n), -1)), perm([[2,1];3:n])]

    for N in 6:20
        S = SN(N)
        S = [S; [perm(randperm(N)) for _ in 1:3]]
        G = PrmGroup(S)
        @test order(G) == factorial(N)
    end

    cube4 = perm.([
       [1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16],
       [1,2,3,4,9,10,11,12,5,6,7,8,13,14,15,16],
       [1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16],
       [16,8,14,6,12,4,10,2,15,7,13,5,11,3,9,1],
       [3,11,1,9,7,15,5,13,4,12,2,10,8,16,6,14]]
       )
    G = PrmGroup(cube4)
    @time StabilizerChain(G)
    @test order(G) == 384

    rubik_gens = [
    perm"( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19)(48)",
    perm"( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35)(48)",
    perm"(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11)(48)",
    perm"(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24)",
    perm"(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27)",
    perm"(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40)"
    ]

    rubik = PrmGroup(rubik_gens)
    @time StabilizerChain(rubik)
    @test order(rubik) == 43252003274489856000 # fits Int128

end
