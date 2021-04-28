function manual_orbit(gens::Vector{GEl}, pt::T, op=^) where {GEl <: GroupElement, T}
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
    @test firstmoved(Perm(5)) == nothing
    @test fixedpoints(perm"(1,2,3,4,5)") isa Vector{Int}
    @test nfixedpoints(perm"(1,2,3,4,5)") == 0
    @test fixedpoints(Perm(Int8[1,2,3,4,5])) isa Vector{Int8}
    @test fixedpoints(perm"(1,3,5)") == [2,4]
    @test nfixedpoints(perm"(1,3,5)") == 2
    @test nfixedpoints(perm"(1,3,5)(6)") == 3

    g = Perm(Int8[1,3,4,2]) # (1)(2,3,4)
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
    @test isone(Perm([1,2,3]))

    @test degree(perm"(1,3)") == 3
    @test degree(perm"(1,3)(5)") == 5

    @test PermutationGroups.emb(perm"(1,3)", 5) != perm"(1,3)"
    @test PermutationGroups.emb(perm"(1,3)", 5) == perm"(1,3)(5)"
    @test PermutationGroups.emb(perm"(1,3)", 5) != perm"(1,3)(6)"
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
    G = SymmetricGroup(SIZE);
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
