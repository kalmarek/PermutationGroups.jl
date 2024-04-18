using Test
using GroupsCore
using PermutationGroups
using BenchmarkTools

function test_perf(G)
    s = 0
    for g in G
        s += 1^g
    end
    return s
end

const GENERATORS = Dict(
    "Sym8_7t" => parse.(Perm{UInt16}, ["($(i), $(i+1))" for i in 1:7]),
    "Sym8_2rp" => [perm"(1,5,6,2,4,8)", perm"(1,3,6)(2,5,7,4)(8)"],
    "cube222" =>
        Perm.([
            [1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16],
            [1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16],
            [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 11, 12, 15, 16],
            [16, 8, 14, 6, 12, 4, 10, 2, 15, 7, 13, 5, 11, 3, 9, 1],
            [3, 11, 1, 9, 7, 15, 5, 13, 4, 12, 2, 10, 8, 16, 6, 14],
        ]),
    "cube333" => [
        perm"( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19)",
        perm"( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35)",
        perm"(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11)",
        perm"(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24)",
        perm"(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27)",
        perm"(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40)",
    ],
    "SL(4,7)" => [
        perm"""
    (  2,  3,  5,  9, 17, 32)(  4,  7, 13, 25, 48, 89)(  6, 11, 21, 40, 76,137)
    (  8, 15, 28, 54, 99,170)( 10, 19, 36, 69,127,113)( 12, 23, 44, 83,148,195)
    ( 14, 27, 52, 96,166,250)( 16, 30, 58,106,179,266)( 18, 34, 65,120,200,290)
    ( 20, 38, 72,131, 39, 74)( 22, 42, 79,142,225,125)( 24, 46, 86,153,220,115)
    ( 26, 50, 92,160,242,232)( 29, 56, 37, 71,130,211)( 31, 60,110,172,260,336)
    ( 33, 63,116,196,287,358)( 35, 67,124,119, 43, 81)( 41, 78,140,223,308, 84)
    ( 45, 80,144,159,117,177)( 47, 53, 93,162,244,194)( 49, 91,158,241,321,264)
    ( 51, 94, 64,118, 68,123)( 55,101,152,121,201,265)( 57,104,178,217,304,171)
    ( 59,108,185,245,324,374)( 61,112,192,285,357,379)( 62,114,191,283,279,351)
    ( 66,122,202, 73, 88,141)( 70,129,209,297,126, 97)( 75,135,216,291,342,247)
    ( 77,139,221,306,334,154)( 82,138,173,262,330,269)( 85,151,235,317,284,205)
    ( 87, 98,168,254,197,289)( 90,157,239,320,259,253)( 95,155,207,294,335,318)
    (100,163,169,256,331,181)(102,174,136,218,204,226)(103,176,175,263,325,156)
    (105,180,268,219,150,228)(107,183,273,257,165,248)(109,187,278,189,261,338)
    (111,190,282,300,366,362)(128,208,296,364,212,145)(132,143,161,243,315,372)
    (133,214,301,288,314,199)(134,206,292,361,376,370)(146,229,313,237,303,369)
    (147,231,164,246,193,286)(149,234,316,373,377,356)(167,252,327,328,305,323)
    (182,271,344,332,238,319)(184,274,186,276,341,382)(198,203,210,222,240,230)
    (213,299,295,215,267,298)(224,309,360,394,375,258)(227,270,343,293,340,307)
    (233,310,255,329,359,339)(236,251,322,311,371,352)(249,326,345,365,386,333)
    (272,346,383,378,275,348)(280,353,389,355,392,350)(302,367,312,337,380,363)
    (347,385,349,387,396,384)(354,391,398)(388,393,390)""",
        perm"""
    (  1,  2,  4,  8, 16, 31, 61,113,193,118,198,183,274, 69,128, 94,164,247,317,309,370,319,367,200,291, 67,125,206,293,362, 63,117, 74,134,215,302,368,392,391,  3,  6, 12, 24, 47, 88,155,238,276,287,359, 99,171,259,335,260,337,381,353,390,  9, 18, 35, 68,126,207,295,363, 40, 77, 54,100,172,261,127,104,179,267,341,283,356,153,236,244,323,218,160,214,262,209,168,255,330,377,361,374,346,384,380, 48, 90, 83,149, 96, 91,159,229,308,216,303,180,269,342,130, 50, 93,163,245,192,116,197, 58,107,184, 76,138,220,305,289,144,157,240, 60,111,191,284, 23, 45, 85,152,176,265, 92,161,110,189,281,355,393, 17, 33, 64,119,199,254,313,364,301,252,328,250,246,325,375,372,344,366,290,360,148,233,315,273,187,279,352, 72,132,213,300, 36, 70, 86,154,237,318,324,326, 65,121, 81,146,230,108,186,277,350,388,397,399,400)
    (  5, 10, 20, 39, 75,136,162,135,217,223,158,129,210,298,282, 11, 22, 43, 82,147,232,292,336,379, 21, 41, 28, 55,102,175,264,286,202,101,173,122,203,271,345,358,327,124,205, 52, 97,167,253,263,339,142,226,311,371,394,211,234,243,185,275,349,190, 89,156,170,258,334,294,343,383,396,357,  7, 14, 15, 29, 57,105,181,270,338,120,174,106,182,272,347,386, 13, 26, 51, 95,165,249, 19, 37, 46, 87,151,139,222,307,348,385,278, 25, 49, 44, 84,150,225,310,321,329,376,248,285,351,369,131,212,178,166,251,201, 42, 80,145,228,297,316,373,320,306,241,221,235,239,242,322,296,140, 78,141,224,304,331,299,365,114,194, 30, 59,109,188,280,354, 32, 62,115,195,268,231,314, 79,143,227,312,196,288, 38, 73,133,208, 71, 56,103,177, 27, 53, 98,169,257,333,137,219,266,340,112, 34, 66,123,204,256,332,378,387,382,395,389,398)""",
    ],
    "direct_product" => [
        inv(Perm(collect([2:401; 1]))),
        Perm([mod(20i, 401) for i in 1:400]),
        perm"""(402,403,404,405,406)""",
        perm"""(402,403)""",
    ],
)

@testset "GAP Docs examples" begin
    @testset "Sym(8) iteration" begin
        G = PermGroup(GENERATORS["Sym8_7t"])
        Ktr = PermGroup(Transversal, GENERATORS["Sym8_2rp"])
        Kstr = PermGroup(SchreierTransversal, GENERATORS["Sym8_2rp"])

        @test factorial(8) ==
              order(Int, G) ==
              order(Int, Ktr) ==
              order(Int, Kstr)

        @test 181440 == test_perf(G) == test_perf(Ktr) == test_perf(Kstr)
    end

    @testset "Rubik cube groups" begin
        RC2tr = PermGroup(Transversal, GENERATORS["cube222"])
        @test order(RC2tr) == 384
        RC2str = PermGroup(SchreierTransversal, GENERATORS["cube222"])
        @test order(RC2str) == 384

        RC3tr = PermGroup(Transversal, GENERATORS["cube333"])
        @test order(RC3tr) == 43252003274489856000
        @test order(RC3tr) == order(Int128, RC3tr) # fits Int128
        RC3str = PermGroup(SchreierTransversal, GENERATORS["cube333"])
        @test order(Int128, RC3str) == 43252003274489856000
    end

    @testset "SL(4,7)" begin
        SL_4_7tr = PermGroup(Transversal, GENERATORS["SL(4,7)"])
        @test order(Int64, SL_4_7tr) == 2317591180800
        SL_4_7str = PermGroup(SchreierTransversal, GENERATORS["SL(4,7)"])
        @test order(Int64, SL_4_7str) == 2317591180800
    end

    @testset "DirectProduct example" begin
        Gtr = PermGroup(Transversal, GENERATORS["direct_product"])
        @test order(Int, Gtr) == 192480
        Gstr = PermGroup(SchreierTransversal, GENERATORS["direct_product"])
        @test order(Int, Gstr) == 192480
    end
end

if !(haskey(ENV, "CI"))
    using BenchmarkTools

    @info "Benchmarking Iteration and Schreier-Sims algorithm"

    begin # iteration Sym(8)
        @info "Iteration over S8 PermGroup (7 transpositions)"
        G = PermGroup(GENERATORS["Sym8_7t"])
        @btime test_perf($G)

        @info "Iteration over K (2 gens) ≅ S8 PermGroup:"
        Ktr = PermGroup(Transversal, GENERATORS["Sym8_2rp"])
        Kstr = PermGroup(SchreierTransversal, GENERATORS["Sym8_2rp"])
        @btime test_perf($Ktr)
        @btime test_perf($Kstr)
    end

    begin # cube groups
        @info "Rubik cube 2×2×2 group:"
        cube222 = GENERATORS["cube222"]
        @btime order(Int, G) setup = (G = PermGroup(Transversal, $cube222)) evals =
            1
        @btime order(Int, G) setup =
            (G = PermGroup(SchreierTransversal, $cube222)) evals = 1

        @info "Schreier-Sims for Rubik cube 3×3×3 group:"
        cube333 = GENERATORS["cube333"]
        @btime order(Int128, G) setup = (G = PermGroup(Transversal, $cube333)) evals =
            1
        @btime order(Int128, G) setup =
            (G = PermGroup(SchreierTransversal, $cube333)) evals = 1
    end

    begin # GAPDocs examples
        @info "Schreier-Sims for SL(4,7):"
        SL_4_7_gens = GENERATORS["SL(4,7)"]
        @btime order(Int64, G) setup =
            (G = PermGroup(Transversal, $SL_4_7_gens)) evals = 1
        @btime order(Int64, G) setup =
            (G = PermGroup(SchreierTransversal, $SL_4_7_gens)) evals = 1

        @info "Schreier-Sims for a direct-product group:"
        S = GENERATORS["direct_product"]
        @btime order(Int, G) setup = (G = PermGroup(Transversal, $S)) evals = 1
        @btime order(Int, G) setup = (G = PermGroup(SchreierTransversal, $S)) evals =
            1

        @info "Iteration over direct-product group of order 192480"
        Gtr = PermGroup(Transversal, S)
        Gstr = PermGroup(SchreierTransversal, S)
        @btime test_perf($Gtr)
        @btime test_perf($Gstr)
    end
end
@testset "SL(4,7)" begin
    SL_4_7tr = PermGroup(Transversal, GENERATORS["SL(4,7)"])
    @test order(Int64, SL_4_7tr) == 2317591180800
    SL_4_7str = PermGroup(SchreierTransversal, GENERATORS["SL(4,7)"])
    @test order(Int64, SL_4_7str) == 2317591180800
end

#=
[ Info: Iteration over S8 PermGroup (7 transpositions)
  3.138 ms (120949 allocations: 6.15 MiB)
[ Info: Iteration over K ≅ S8 PermGroup:
  0.000039 seconds (562 allocations: 32.742 KiB)
  3.121 ms (120949 allocations: 6.15 MiB)
  0.000075 seconds (1.16 k allocations: 70.273 KiB)
  5.809 ms (201815 allocations: 11.09 MiB)
[ Info: Testing and benchmarking Schreier-Sims algorithm
[ Info: Rubik cube 2×2×2 group:
  18.565 μs (442 allocations: 29.43 KiB)
  52.597 μs (1148 allocations: 79.07 KiB)
[ Info: Schreier-Sims for Rubik cube 3×3×3 group:
  1.374 ms (24680 allocations: 2.31 MiB)
  6.975 ms (111829 allocations: 10.81 MiB)
[ Info: Schreier-Sims for SL(4,7):
  4.393 ms (25970 allocations: 10.22 MiB)
  53.474 ms (189388 allocations: 76.13 MiB)
[ Info: Schreier-Sims for a direct-product group:
  1.228 ms (9104 allocations: 3.15 MiB)
  58.843 ms (219502 allocations: 77.17 MiB)
[ Info: Iteration over direct-product group of order 192480
  0.061710 seconds (585.73 k allocations: 184.881 MiB, 9.19% gc time)
  0.201417 seconds (1.20 M allocations: 441.331 MiB, 6.48% gc time)
=#
