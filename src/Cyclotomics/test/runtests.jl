using Test
using Cyclotomics
import Cyclotomics.Cyclotomic

@testset "Cyclotomics" begin

    @testset "zumbroich" begin

        @test Cyclotomics._forbidden_residues(45, 5, 1) == BitSet([ 0 ])
        @test Cyclotomics._forbidden_residues(45, 3, 2) == BitSet([ 4, 0, 5 ])

        @test last(first(Cyclotomics.ForbiddenResidues(9))) isa BitSet
        fb = Cyclotomics.ForbiddenResidues(9)
        @test 0 in fb
        @test 1 in fb
        @test !(2 in fb)

        @test Cyclotomics.zumbroich(9) == [ 2, 3, 4, 5, 6, 7 ]

        @test !any(in(fb), Cyclotomics.zumbroich(9))

        @test Cyclotomics.zumbroich_plain(8) ==
            Cyclotomics.zumbroich(8) == [ 0, 1, 2, 3 ]
        @test Cyclotomics.zumbroich_plain(9) ==
            Cyclotomics.zumbroich(9) == [ 2, 3, 4, 5, 6, 7 ]

        @test Cyclotomics.zumbroich(45) == [ 1, 2, 3, 6, 7, 8, 11, 12, 16, 17, 19, 21, 24, 26, 28, 29, 33, 34, 37, 38, 39, 42, 43, 44 ]

        @test all(Cyclotomics.zumbroich_plain(i) ==
            Cyclotomics.zumbroich(i) ==
            Cyclotomics.zumbroich_direct(i) for i in 1:5000)
    end

    @testset "elementary ops" begin
        @test E(5) isa Number
        @test E(5) isa Cyclotomic{Int}
        @test E(5, 0) isa Cyclotomic{Int}
        @test E(5, 6) isa Cyclotomic{Int}

        E₅ = E(5)
        @test valtype(E₅) == Int
        @test Cyclotomic{Float64}(E₅) isa Cyclotomic{Float64}

        E₅fl = Cyclotomic{Float64}(E₅)
        @test valtype(E₅fl) == Float64

        @test similar(E₅) isa Cyclotomic{Int}
        @test similar(E₅fl) isa Cyclotomic{Float64}

        @test zero(E₅) isa Cyclotomic{Int}
        @test zero(E₅fl) isa Cyclotomic{Float64}

        @test one(E₅) isa Cyclotomic{Int}
        @test one(E₅fl) isa Cyclotomic{Float64}

        @test deepcopy(E₅).coeffs !== E₅.coeffs
        @test deepcopy(E₅).coeffs == E₅.coeffs
    end

    @testset "io strings" begin
        @test string(E(5)) == " +1*E(5)^1"
        @test string(-E(5)) == " -1*E(5)^1"

        using Base.Meta
        x = E(5) + 2E(5)^2
        @test string(x) == " +1*E(5)^1 +2*E(5)^2"
        @test eval(Base.Meta.parse(string(x))) == x

        x = E(5) - 2E(5)^2
        @test string(x) == " +1*E(5)^1 -2*E(5)^2"
        @test eval(Base.Meta.parse(string(x))) == x

        x = -E(5) + 2E(5)^2
        @test string(x) == " -1*E(5)^1 +2*E(5)^2"
        @test eval(Base.Meta.parse(string(x))) == x
    end

    @testset "indexing" begin
        x = E(5)
        @test x[0] == 0
        @test x[1] == 1
        @test all(iszero, x[2:5])
        @test x[-1] == 0
        @test x[-4] == 1
        @test x[6] == 1

        @test setindex!(x, 1, 2) == 1
        x[2] = 3
        @test x[0] == 0
        @test x[1] == 1
        @test x[2] == 3
        @test x[-3] == 3
        @test x[7] == 3
    end

    @testset "aritmetic: +, -, module: *, //" begin
        x = E(5)
        y = E(5,2)

        @test 2x isa Cyclotomic{Int}
        @test 2.0x isa Cyclotomic{Float64}
        @test div(x, 2) isa Cyclotomic{Int}
        @test x//2 isa Cyclotomic{Rational{Int}}
        @test x/2.0 isa Cyclotomic{Float64}
        @test x/2 isa Cyclotomic{Float64}

        @test x+2y isa Cyclotomic{Int}
        xy = x+2y
        @test xy[0] == 0
        @test xy[1] == 1
        @test xy[2] == 2
        @test all(iszero, xy[3:5])

        @test 2.0xy[0] isa Float64
        @test 2.0xy[1] == 2

        @test (xy - 2y)[1] == 1
        @test all(iszero, (xy - 2y)[2:5])

        @test 1 + x isa Cyclotomic{Int}
        @test (1+x)[0] == 1
        @test x + 1 isa Cyclotomic{Int}
        @test 2.0 + x isa Cyclotomic{Float64}
        @test x + 2.0 isa Cyclotomic{Float64}
        @test (x+2.0)[0] == 2.0
    end

    @testset "*, powering" begin
        x = E(5)
        y = E(5,2)

        w = x*(x+2y)
        @test x*y isa Cyclotomic{Int}
        @test (x*y)[1] == 0
        @test (x*y)[2] == 0
        @test (x*y)[3] == 1
        w = (1+x)*(1+y)
        @test w[0] == 1
        @test w[1] == 1
        @test w[2] == 1
        @test w[3] == 1

        @test (x^2)[1] == 0
        @test (x^2)[2] == 1

        @test ((1+x)^2)[0] == 1
        @test ((1+x)^2)[1] == 2
        @test ((1+x)^2)[2] == 1
        @test ((1+x)^2)[3] == 0

        @test ((1+x^3)^2)[0] == 1
        @test ((1+x^3)^2)[1] == 1
        @test ((1+x^3)^2)[3] == 2
    end

    @testset "normal form" begin
        x = E(45) + E(45, 5);
        x.coeffs
        @test deepcopy(x) !== x

        y = Cyclotomics.normalform(x);

        e = E(45)
        w = e + e^2 + e^8 + e^11 + e^17 + e^26 + e^29 + e^38 + e^44

        @test w == y
        @test x.coeffs != y.coeffs
        @test y.coeffs == w.coeffs

        @test deepcopy(x) == deepcopy(y)
        @test x !== deepcopy(x)
        @test x.coeffs === copy(x).coeffs

        @test hash(deepcopy(x)) == hash(deepcopy(y))
        @test length(Set([deepcopy(x), deepcopy(y), deepcopy(x), deepcopy(y)])) == 1

        @test iszero(1 + x - x - 1)

        @test isreal(1+x-x)
        @test isone(sum(-E(5)^i for i in 1:4))
    end

    @testset "embedding" begin
        let x = E(45)^5 + E(45)^10
            @test conductor(Cyclotomics.reduced_embedding(x)) == 9
            y = Cyclotomics.reduced_embedding(x)
            @test y == E(9)^2 - E(9)^4 - E(9)^7
            Cyclotomics.normalform!(x)
            @test_broken conductor(Cyclotomics.reduced_embedding(x)) == 9
        end
    end

    @testset "tests against GAP" begin
        @test E(9) == -E(9)^4 - E(9)^7
        @test E(9)^3 == E(3)
        @test E(6) == -E(3)^2

        @test E(45)^ 4 == -E(45)^19-E(45)^34
        @test E(45)^13 == -E(45)^28-E(45)^43
        @test E(45)^14 == -E(45)^29-E(45)^44
        @test E(45)^22 == -E(45)^ 7-E(45)^37

        @test_broken E(5) + E(3) == -E(15)^2-2*E(15)^8-E(15)^11-E(15)^13-E(15)^14
        @test (E(5) + E(5)^4) ^ 2 == -2*E(5)-E(5)^2-E(5)^3-2*E(5)^4
        @test_broken E(5) / E(3) == E(15)^13
        @test_broken E(5) * E(3) == E(15)^8
    end

end
