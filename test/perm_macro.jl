@testset "perm_str macro" begin
    PG = PermutationGroups
    @test PG.Perms._parse_cycles("()") == Int[]
    @test PG.Perms._parse_cycles("(1)(2)(3)") == [[1], [2], [3]]
    @test PG.Perms._parse_cycles("(1)(2,3)") == [[1], [2, 3]]
    @test PG.Perms._parse_cycles("(1)(\n2, 3)") == [[1], [2, 3]]
    @test PG.Perms._parse_cycles("(3,2,1)(4,5)") == [[3, 2, 1], [4, 5]]
    @test_throws ArgumentError PG.Perms._parse_cycles("(a,b)")
    @test_throws ArgumentError PG.Perms._parse_cycles("(1 2)")

    s = """
 ( 1, 22,73,64,78,81,  24 ,89,90,54,51,82,91,53, 18
  ,38,19,52,44,77,62,95,94,50,43,42,
 10,67,87,60,36,12)(2,57,34,88)(3,92,76,17,99,96,30,55,45,41,98)(4,56,59,97,49,
 21,15,9,26,86,83,29,27,66,6,58,28,5,68,40,72,7,84,93,39,79,23,46,63,32,61,100,
 11)(8,80,71,75,35,14,85,25,20,70,65,16,48,47,37,74,33,13,31,69)
 """

    s2 = """
  (1,22,73,64,78,81,24,89,90,54,51,82,91,53,18,38,19,52,44,77,62,95,94,50,43,42,\n10,67,87,60,36,12)(2,57,34,88)(3,92,76,17,99,96,30,55,45,41,98)(4,56,59,97,49,\n21,15,9,26,86,83,29,27,66,6,58,28,5,68,40,72,7,84,93,39,79,23,46,63,32,61,100,\n11)(8,80,71,75,35,14,85,25,20,70,65,16,48,47,37,74,33,13,31,69)
  """
    @test PG.Perms._parse_cycles(s) == PG.Perms._parse_cycles(s2)
    P = Perm{UInt8}
    images = UInt8[0x02, 0x03, 0x01]
    @test Meta.parse(Perm{UInt8}, "(1,2,3)(5)(10)") == :($P($images))
    @test Meta.parse(Perm{UInt32}, "(1,2,3)(5)(10)") ==
          :($(Perm{UInt32})($(UInt32.(images))))

    @test perm"(1,2,3)(5)(10)" isa PG.AbstractPermutation
    @test perm"(1,2,3)(5)(10)" isa PG.Perm
    @test perm"(1,2,3)(5)(10)" isa PG.Perm{UInt16}
    @test parent(perm"(1,2,3)(5)(10)") == PG.Perms.InfinitePermGroup()
    @test degree(perm"(1,2,3)(5)(10)") == 3
    @test perm"(1,2,3,4,5)" == Perm([2, 3, 4, 5, 1])
    @test perm"(3,2,1)(4,5)" == Perm([3, 1, 2, 5, 4])
end
