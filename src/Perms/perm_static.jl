function __ntuple_u8(itr, ::Val{N}, d = length(itr)) where {N}
    @inbounds t =
        ntuple(i -> i ≤ d ? UInt8(itr[i]) : UInt8(i), N)::NTuple{N,UInt8}
    return t
end

function __degree(images::Tuple)
    n = length(images)
    for i in n:-1:1
        images[i] ≠ i && return i
    end
    return 0
end

struct SPerm{N} <: AP.AbstractPermutation
    images::NTuple{N,UInt8}
    degree::UInt8
end

function SPerm{N}(
    images::AbstractVector{<:Integer};
    check::Bool = true,
) where {N}
    @assert length(images) ≤ N "vector too large for SPerm"
    @assert N ≤ typemax(UInt8)
    if check && !isperm(images)
        throw(ArgumentError("Provided images are not permutation!"))
    end
    deg = __degree(images)
    dimages = @view images[firstindex(images):deg]

    return SPerm{N}(__ntuple_u8(dimages, Val{N}(), deg), deg)
end

function Base.:^(n::Integer, σ::SPerm)
    return 1 ≤ n ≤ AP.degree(σ) ? oftype(n, @inbounds σ.images[n]) : n
end
AP.degree(σ::SPerm) = σ.degree
AP.inttype(::Type{<:SPerm}) = UInt8
AP.__unsafe_image(n::Integer, σ::SPerm) = oftype(n, @inbounds σ.images[n])

function Base.inv(σ::SPerm{N}) where {N}
    return SPerm{N}(invperm(σ.images), AP.degree(σ))
end

function Base.:*(σ::SPerm{N}, τ::SPerm{M}) where {N,M}
    K = max(N, M)
    t = ntuple(i -> (i^σ)^τ, K)
    return SPerm{K}(t, __degree(t))
end

macro sperm8_str(str)
    return :(Base.parse($(esc(SPerm{8})), $str))
end
