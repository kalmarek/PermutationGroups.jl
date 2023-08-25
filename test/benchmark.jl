using Test
using GroupsCore
using PermutationGroups
using BenchmarkTools

const BENCHMARK_TIME = !(haskey(ENV, "CI"))

function test_perf(G)
    s = 0
    for g in G
        s += 1^g
    end
    return s
end

@testset "GAP Docs examples$(BENCHMARK_TIME ? "/benchmark" : "")" begin
    @testset "test_perf$(BENCHMARK_TIME ? "/benchmark" : "") iteration" begin
        import AbstractAlgebra
        SymmetricGroup = AbstractAlgebra.SymmetricGroup

        Base.:^(i::Integer, p::AbstractAlgebra.Perm) = oftype(i, p.d[i])

        G = SymmetricGroup(8)
        K = PermGroup(
            Transversal,
            [perm"(1,5,6,2,4,8)", perm"(1,3,6)(2,5,7,4)(8)"],
        )
        @test order(Int, K) == order(Int, G)
        @test test_perf(G) == test_perf(K) == 181440
        K = PermGroup(
            SchreierTransversal,
            [perm"(1,5,6,2,4,8)", perm"(1,3,6)(2,5,7,4)(8)"],
        )
        @test order(Int, K) == order(Int, G)
        @test test_perf(G) == test_perf(K) == 181440

        if BENCHMARK_TIME
            @info "Native iteration over S8 group:"
            @btime test_perf($G)
            # 1.303 ms (80644 allocations: 6.77 MiB)
            @info "Iteration over K ≅ S8 PermGroup:"
            K = PermGroup(
                Transversal,
                [perm"(1,5,6,2,4,8)", perm"(1,3,6)(2,5,7,4)(8)"],
            )
            # 3.285 ms (125987 allocations: 6.46 MiB)
            @time order(Int, K)
            @btime test_perf($K)
            K = PermGroup(
                SchreierTransversal,
                [perm"(1,5,6,2,4,8)", perm"(1,3,6)(2,5,7,4)(8)"],
            )
            @time order(Int, K)
            @btime test_perf($K)
            # 6.833 ms (206127 allocations: 11.35 MiB)
        end
    end

    if BENCHMARK_TIME
        @info "Testing and benchmarking Schreier-Sims algorithm"
    end
    @testset "Rubik cube groups" begin
        cube222 =
            Perm.([
                [1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16],
                [1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 13, 14, 15, 16],
                [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 11, 12, 15, 16],
                [16, 8, 14, 6, 12, 4, 10, 2, 15, 7, 13, 5, 11, 3, 9, 1],
                [3, 11, 1, 9, 7, 15, 5, 13, 4, 12, 2, 10, 8, 16, 6, 14],
            ])
        G = PermGroup(Transversal, cube222)
        @test order(G) == 384
        G = PermGroup(SchreierTransversal, cube222)
        @test order(G) == 384
        if BENCHMARK_TIME
            @info "Rubik cube 2×2×2 group:"
            @btime order(Int, G) setup = (G = PermGroup(Transversal, $cube222)) evals =
                1
            @btime order(Int, G) setup =
                (G = PermGroup(SchreierTransversal, $cube222)) evals = 1
            # 24.767 μs (533 allocations: 44.37 KiB)
            # 63.420 μs (1239 allocations: 109.91 KiB)
        end

        cube333 = [
            perm"( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19)",
            perm"( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35)",
            perm"(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11)",
            perm"(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24)",
            perm"(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27)",
            perm"(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40)",
        ]

        G = PermGroup(Transversal, cube333)
        @test order(G) == 43252003274489856000
        @test order(G) == order(Int128, G) # fits Int128
        G = PermGroup(SchreierTransversal, cube333)
        @test order(Int128, G) == 43252003274489856000
        if BENCHMARK_TIME
            @info "Schreier-Sims for Rubik cube 3×3×3 group:"
            @btime order(Int128, G) setup =
                (G = PermGroup(Transversal, $cube333)) evals = 1
            @btime order(Int128, G) setup =
                (G = PermGroup(SchreierTransversal, $cube333)) evals = 1
            # 1.975 ms (26628 allocations: 2.59 MiB)
            # 8.836 ms (116353 allocations: 11.42 MiB)
        end
    end

    @testset "SL(4,7)" begin
        a = perm"""
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
    (347,385,349,387,396,384)(354,391,398)(388,393,390)"""

        b =
            perm"""(  1,  2,  4,  8, 16, 31, 61,113,193,118,198,183,274, 69,128, 94,164,247,317,309,370,319,367,200,291, 67,125,206,293,362, 63,117, 74,134,215,302,368,392,391,  3,  6, 12, 24, 47, 88,155,238,276,287,359, 99,171,259,335,260,337,381,353,390,  9, 18, 35, 68,126,207,295,363, 40, 77, 54,100,172,261,127,104,179,267,341,283,356,153,236,244,323,218,160,214,262,209,168,255,330,377,361,374,346,384,380, 48, 90, 83,149, 96, 91,159,229,308,216,303,180,269,342,130, 50, 93,163,245,192,116,197, 58,107,184, 76,138,220,305,289,144,157,240, 60,111,191,284, 23, 45, 85,152,176,265, 92,161,110,189,281,355,393, 17, 33, 64,119,199,254,313,364,301,252,328,250,246,325,375,372,344,366,290,360,148,233,315,273,187,279,352, 72,132,213,300, 36, 70, 86,154,237,318,324,326, 65,121, 81,146,230,108,186,277,350,388,397,399,400)
        (  5, 10, 20, 39, 75,136,162,135,217,223,158,129,210,298,282, 11, 22, 43, 82,147,232,292,336,379, 21, 41, 28, 55,102,175,264,286,202,101,173,122,203,271,345,358,327,124,205, 52, 97,167,253,263,339,142,226,311,371,394,211,234,243,185,275,349,190, 89,156,170,258,334,294,343,383,396,357,  7, 14, 15, 29, 57,105,181,270,338,120,174,106,182,272,347,386, 13, 26, 51, 95,165,249, 19, 37, 46, 87,151,139,222,307,348,385,278, 25, 49, 44, 84,150,225,310,321,329,376,248,285,351,369,131,212,178,166,251,201, 42, 80,145,228,297,316,373,320,306,241,221,235,239,242,322,296,140, 78,141,224,304,331,299,365,114,194, 30, 59,109,188,280,354, 32, 62,115,195,268,231,314, 79,143,227,312,196,288, 38, 73,133,208, 71, 56,103,177, 27, 53, 98,169,257,333,137,219,266,340,112, 34, 66,123,204,256,332,378,387,382,395,389,398)"""

        SL_4_7 = PermGroup(Transversal, [a, b])
        @test order(Int64, SL_4_7) == 2317591180800
        SL_4_7 = PermGroup(SchreierTransversal, [a, b])
        @test order(Int64, SL_4_7) == 2317591180800
        if BENCHMARK_TIME
            @info "Schreier-Sims for SL(4,7):"
            @btime order(Int64, G) setup =
                (G = PermGroup(Transversal, $([a, b]))) evals = 1
            @btime order(Int64, G) setup =
                (G = PermGroup(SchreierTransversal, $([a, b]))) evals = 1
            # gap> G := Group([a,b]);; StabChain(G);; time;
            # ~35ms
            # 10.258 ms (30072 allocations: 11.49 MiB)
            # 98.011 ms (240549 allocations: 102.46 MiB)
        end
    end

    @testset "DirectProduct example" begin
        a = perm"""
        (  1,401,400,399,398,397,396,395,394,393,392,391,390,389,388,387,386,
         385,384,383,382,381,380,379,378,377,376,375,374,373,372,371,370,369,
         368,367,366,365,364,363,362,361,360,359,358,357,356,355,354,353,352,
         351,350,349,348,347,346,345,344,343,342,341,340,339,338,337,336,335,
         334,333,332,331,330,329,328,327,326,325,324,323,322,321,320,319,318,
         317,316,315,314,313,312,311,310,309,308,307,306,305,304,303,302,301,
         300,299,298,297,296,295,294,293,292,291,290,289,288,287,286,285,284,
         283,282,281,280,279,278,277,276,275,274,273,272,271,270,269,268,267,
         266,265,264,263,262,261,260,259,258,257,256,255,254,253,252,251,250,
         249,248,247,246,245,244,243,242,241,240,239,238,237,236,235,234,233,
         232,231,230,229,228,227,226,225,224,223,222,221,220,219,218,217,216,
         215,214,213,212,211,210,209,208,207,206,205,204,203,202,201,200,199,
         198,197,196,195,194,193,192,191,190,189,188,187,186,185,184,183,182,
         181,180,179,178,177,176,175,174,173,172,171,170,169,168,167,166,165,
         164,163,162,161,160,159,158,157,156,155,154,153,152,151,150,149,148,
         147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,
         130,129,128,127,126,125,124,123,122,121,120,119,118,117,116,115,114,
         113,112,111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97,
          96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80,
          79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63,
          62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46,
          45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29,
          28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12,
          11, 10,  9,  8,  7,  6,  5,  4,  3,  2)"""
        b = perm"""
        (  1, 20,400,381)(  2, 40,399,361)(  3, 60,398,341)(  4, 80,397,321)
        (  5,100,396,301)(  6,120,395,281)(  7,140,394,261)(  8,160,393,241)
        (  9,180,392,221)( 10,200,391,201)( 11,220,390,181)( 12,240,389,161)
        ( 13,260,388,141)( 14,280,387,121)( 15,300,386,101)( 16,320,385, 81)
        ( 17,340,384, 61)( 18,360,383, 41)( 19,380,382, 21)( 22, 39,379,362)
        ( 23, 59,378,342)( 24, 79,377,322)( 25, 99,376,302)( 26,119,375,282)
        ( 27,139,374,262)( 28,159,373,242)( 29,179,372,222)( 30,199,371,202)
        ( 31,219,370,182)( 32,239,369,162)( 33,259,368,142)( 34,279,367,122)
        ( 35,299,366,102)( 36,319,365, 82)( 37,339,364, 62)( 38,359,363, 42)
        ( 43, 58,358,343)( 44, 78,357,323)( 45, 98,356,303)( 46,118,355,283)
        ( 47,138,354,263)( 48,158,353,243)( 49,178,352,223)( 50,198,351,203)
        ( 51,218,350,183)( 52,238,349,163)( 53,258,348,143)( 54,278,347,123)
        ( 55,298,346,103)( 56,318,345, 83)( 57,338,344, 63)( 64, 77,337,324)
        ( 65, 97,336,304)( 66,117,335,284)( 67,137,334,264)( 68,157,333,244)
        ( 69,177,332,224)( 70,197,331,204)( 71,217,330,184)( 72,237,329,164)
        ( 73,257,328,144)( 74,277,327,124)( 75,297,326,104)( 76,317,325, 84)
        ( 85, 96,316,305)( 86,116,315,285)( 87,136,314,265)( 88,156,313,245)
        ( 89,176,312,225)( 90,196,311,205)( 91,216,310,185)( 92,236,309,165)
        ( 93,256,308,145)( 94,276,307,125)( 95,296,306,105)(106,115,295,286)
        (107,135,294,266)(108,155,293,246)(109,175,292,226)(110,195,291,206)
        (111,215,290,186)(112,235,289,166)(113,255,288,146)(114,275,287,126)
        (127,134,274,267)(128,154,273,247)(129,174,272,227)(130,194,271,207)
        (131,214,270,187)(132,234,269,167)(133,254,268,147)(148,153,253,248)
        (149,173,252,228)(150,193,251,208)(151,213,250,188)(152,233,249,168)
        (169,172,232,229)(170,192,231,209)(171,212,230,189)(190,191,211,210)"""
        c = perm"""(402,403,404,405,406)"""
        d = perm"""(402,403)"""

        G = PermGroup(Transversal, [a, b, c, d])
        @test order(Int, G) == 192480
        G = PermGroup(SchreierTransversal, [a, b, c, d])
        @test order(Int, G) == 192480
        if BENCHMARK_TIME
            @info "Schreier-Sims for a direct-product group:"
            @btime order(Int, G) setup =
                (G = PermGroup(Transversal, $([a, b, c, d]))) evals = 1
            @btime order(Int, G) setup =
                (G = PermGroup(SchreierTransversal, $([a, b, c, d]))) evals = 1
            # gap> G := Group([a,b,c,d]);; StabChain(G);; time;
            # ~35ms
            # 2.759 ms (9133 allocations: 3.18 MiB)
            # 83.738 ms (221318 allocations: 77.73 MiB)
            G = PermGroup(Transversal, [a, b, c, d])
            @time test_perf(G)
            H = PermGroup(SchreierTransversal, [a, b, c, d])
            @time test_perf(H)
            # 0.146197 seconds (622.65 k allocations: 216.433 MiB, 3.77% gc time)
            # 0.388392 seconds (1.24 M allocations: 473.432 MiB, 3.57% gc time)
        end
    end
end
