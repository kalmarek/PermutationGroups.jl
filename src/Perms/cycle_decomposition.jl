struct CycleDecomposition{T<:Integer}
    cycles::Vector{T} # cycles, concatenated
    cycles_ptrs::Vector{T} # pointers to the starts of the cycles
end

Base.length(cd::CycleDecomposition) = length(cd.cycles_ptrs) - 1
function Base.eltype(::Type{CycleDecomposition{T}}) where {T}
    return SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}
end

function Base.iterate(cd::CycleDecomposition, state = 1)
    state == length(cd.cycles_ptrs) && return nothing
    from = cd.cycles_ptrs[state]
    to = cd.cycles_ptrs[state+1] - 1
    return @view(cd.cycles[from:to]), state + 1
end

function Base.show(io::IO, cd::CycleDecomposition)
    print(io, "Cycle Decomposition: ")
    for c in cd
        print(io, '(')
        join(io, c, ',')
        print(io, ')')
    end
end
