function _parse_cycles(str::AbstractString)
    cycles = Vector{Vector{Int}}()
    if occursin(r"\d\s+\d", str)
        throw(ArgumentError("parse string as cycles: spaces between digits"))
    end
    str = replace(str, r"\s+" => "")
    str = replace(str, "()" => "")
    cycle_regex = r"\(\d+(,\d+)*\)?"
    parsed_size = 0
    for m in eachmatch(cycle_regex, str)
        cycle_str = m.match
        parsed_size += sizeof(cycle_str)
        cycle = [parse(Int, a) for a in split(cycle_str[2:end-1], ",")]
        push!(cycles, cycle)
    end
    if parsed_size != sizeof(str)
        throw(
            ArgumentError(
                "parse string as cycles: parsed size differs from string",
            ),
        )
    end
    return cycles
end

function Meta.parse(
    ::Type{P},
    str::AbstractString,
) where {P<:AbstractPermutation}
    cycles = _parse_cycles(str)
    deg = mapreduce(maximum, max, cycles; init = 1)
    images = Vector{inttype(P)}(undef, deg)
    for idx in 1:deg
        k = idx
        for cycle in cycles
            i = findfirst(==(k), cycle)
            k = isnothing(i) ? k : cycle[mod1(i + 1, length(cycle))]
        end
        images[idx] = k
    end
    return :($P($images))
end
"""
    perm"..."

String macro to parse cycles into `Perm{UInt16}`.

Strings for the output of GAP could be copied directly into `perm"..."`.
Cycles of length `1` are not necessary, but can be included.

# Examples:
```jldoctest
julia> p = perm"(1,3)(2,4)"
(1,3)(2,4)

julia> typeof(p)
Perm{Int64}

julia> degree(p)
4

julia> q = perm"(1,3)(2,4)(10)"
(1,3)(2,4)

julia> degree(q)
4
```
"""
macro perm_str(str)
    return Meta.parse(Perm{inttype(Perm)}, str)
end
