_inv_of(cc::AbstractVector) = _inv_of(Int16, cc)

function _inv_of(::Type{I}, cc::AbstractVector) where I<:Integer
    inv_of = zeros(I, size(cc))
    for (i, c) in enumerate(cc)
        g = inv(first(c))
        inv_of[i] = something(findfirst(k -> g in k, cc), 0)
        @assert !iszero(inv_of[i]) "Could not find the conjugacy class of $g."
    end
    return inv_of
end

Base.eltype(::AbstractClassFunction{T}) where T = T

function Base.getindex(χ::Character, i::Integer)
    @boundscheck checkbounds(χ.vals, abs(i))
    if i < 0
        return @inbounds χ.vals[χ.inv_of[abs(i)]]
    else
        return @inbounds χ.vals[i]
    end
end

function (χ::AbstractClassFunction)(g::GroupElem)
    for (i,cc) in enumerate(χ.cc)
        g ∈ cc && return χ[i]
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
