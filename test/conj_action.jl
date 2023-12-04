@testset "Conjugation action" begin
    S = [perm"(1,2)(4)", perm"(1,2,3,4)"]

    pt = perm"(4)"
    @test PG.Orbit(pt, S) isa Orbit
    @test length(Orbit(pt, S)) == 1

    pt = perm"(1,2)(4)"
    @test Orbit(pt, S) isa Orbit
    @test length(Orbit(pt, S)) == 6
    @test all(==(AP.permtype(pt)), AP.permtype.(Orbit(pt, S)))

    pt = perm"(1,2)(4)"
    @test Transversal(pt, S) isa Transversal
    tr = Transversal(pt, S)
    @test all(pt^tr[p] == p for p in tr)

    pt = perm"(2,3,4)"
    @test SchreierTransversal(pt, S) isa SchreierTransversal
    tr = SchreierTransversal(pt, S)
    @test all(pt^tr[p] == p for p in tr)
end
