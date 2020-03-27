Base.eltype(::AbstractCharacter{T}) where T = T

function (χ::Character)(g::GroupElem)
    for (i,cc) in enumerate(χ.cc)
        g ∈ cc && return χ.vals[i]
    end
    throw(DomainError(g, "element does not belong to conjugacy classes of χ"))
end

function LinearAlgebra.dot(v::AbstractCharacter{T}, w::AbstractCharacter{T}) where T<:FinFieldElem
    # TODO: @assert v.cc == w.cc

    F = parent(first(v.vals)) # TODO make something better here
    ord = F(sum(length, v.cc))
    val = mapreduce(
        (i, cc) -> v.vals[i]*w.vals[i]*inv(length(cc)),
        +,
        enumerate(v.cc),
        init = F(0)
    )
    val *= inv(ord)
    return sqrt(val)
end
