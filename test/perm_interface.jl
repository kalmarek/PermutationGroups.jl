module APerms

using Test
import PermutationGroups
const PG = PermutationGroups

export APerm

struct APerm <: PG.AbstractPermutation
    images::Vector{Int}
    APerm(images, check::Bool = true) = new(images) # no checks :)
end

@testset "Implementing AbstractPermutation interface" begin
    @test one(APerm) isa PG.AbstractPermutation

    @test_throws String PG.degree(one(APerm))

    function PG.degree(p::APerm)
        return something(findlast(i -> p.images[i] ≠ i, eachindex(p.images)), 1)
    end

    @test PG.degree(one(APerm)) == 1

    @test_throws String 5^one(APerm)

    function Base.:^(i::Integer, p::APerm)
        return 1 ≤ i ≤ PG.degree(p) ? oftype(i, p.images[i]) : i
    end

    @test 5^one(APerm) == 5

    @test PG.Perms.inttype(one(APerm)) == UInt32
    # but actually it'd be better to have it as Int64
    one(APerm)
    k1 = @allocated one(APerm)
    PG.Perms.inttype(::Type{APerm}) = Int

    one(APerm)
    k2 = @allocated one(APerm)
    @test k2 < k1
end

end
